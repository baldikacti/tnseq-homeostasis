# cython: boundscheck = False
# cython: initializedcheck = False
# cython: wraparound = False
# cython: cdivision = True
# cython: language_level = 3

import sys
import numpy as np
import numpy.random as rn
# import scipy.stats as ss
cimport numpy as np
from libc.math cimport sqrt, log1p, exp, log

from cython.parallel import parallel, prange

from apf.base.apf cimport APF
from apf.base.sample cimport _sample_gamma, _sample_poisson, _sample_crt, _sample_categorical, _sample_logcategorical
from apf.base.cyutils cimport _logsumexp, _logpmf_poisson

cdef extern from "gsl/gsl_rng.h" nogil:
    ctypedef struct gsl_rng:
        pass

def exit_if(func_output, func_desc):
    if func_output is not None:
        if not func_output:
            sys.exit('Error in %s. Exiting.' % func_desc)

cdef class GaP_HDP_Seq(APF):

    cdef:
        int n_loci, n_bases, n_genotypes, n_samples, n_subjects, n_subpops
        double gam, beta, shp_gam, rte_gam, shp_beta, rte_beta
        int[:,:,::1] T_GB
        double[:,:,::1] P_XLG
        double[:,::1] Theta_SK, Theta_IK, A_LG
        double[::1] Theta_K
        int[:,::1] Z_LK
        long[::1] H_K
        long[:,::1] H_IK
        long[::1] n_biopsies_I
        double[:,::1] H_T
        double[::1] phi_lk

    def __init__(self, data_SM, n_loci, list n_biopsies_I, int n_subpops,
                 shp_gam=0.1, rte_gam=1., shp_beta=10., rte_beta=10.,
                 object seed=None, object n_threads=None):

    	############### H -> Z_LK, G -> Theta, K, IK, SK ########################
        # self.T_GB = np.asarray([[0.001], [1000]], dtype='double') #### change
        # self.n_genotypes = self.T_GB.shape[0]
        # self.n_bases = self.T_GB.shape[1]
        self.n_genotypes = 2
        self.n_bases = 1
        self.A_LG = np.repeat([[0.5, 0.5]], n_loci, axis=0) ##### change
        # self.A_LG = np.repeat([[0.5, 0.5]], n_loci, axis=0) ##### change
        self.n_loci = n_loci                                      # L
        self.n_subjects = len(n_biopsies_I)                       # I
        self.n_biopsies_I = np.array(n_biopsies_I).astype(int)    # J_i
        self.n_samples = np.sum(n_biopsies_I)                     # S
        self.n_subpops = n_subpops                                # K/Q

        data_shp = (self.n_samples, self.n_loci * self.n_bases)  # (S, LxB)
        core_shp = (self.n_subpops,)

        # self.T_GB = np.zeros((self.n_loci, self.n_genotypes, self.n_bases), dtype='double')
        with open('threshold.npy', 'rb') as f:
            self.T_GB = np.load(f).reshape(-1,2,1).astype(np.int32) 
        # for l in range(self.n_loci):
        #     for g in range(self.n_genotypes):
        #         if g == 0:
        #             self.T_GB[l, g, 0] = 0.001
        #         elif g == 1:
        #             # self.T_GB[l, g, 0] = np.percentile(data_SM[:, l], 50)
        #             mid = np.max(data_SM[:, l]) / 10
        #             if mid > 100:
        #                 self.T_GB[l, g, 0] = mid
        #             else:
        #                 self.T_GB[l, g, 0] = 100
        #             # self.T_GB[l, g, 0] = np.max(data_SM[:, l])
        #             # self.T_GB[l, g, 0] = np.mean(data_SM[:, l])
        #             # self.T_GB[l, g, 0] = 100
        #         else:
        #             # self.T_GB[l, g, 0] = np.percentile(data_SM[:, l], 75)
        #             count_max = np.max(data_SM[:, l])
        #             if count_max > 1000:
        #                 self.T_GB[l, g, 0] = np.max(data_SM[:, l])
        #             else:
        #                 self.T_GB[l, g, 0] = 1000
        #             # self.T_GB[l, g, 0] = 1000


        super(GaP_HDP_Seq, self).__init__(data_shp=data_shp,
                                          core_shp=core_shp,
                                          eps=0.1,
                                          binary=0,
                                          mtx_is_dirichlet=[1],
                                          seed=seed,
                                          n_threads=n_threads)
        # Params
        self.shp_gam = self.param_list['shp_gam'] = shp_gam
        self.rte_gam = self.param_list['rte_gam'] = rte_gam
        self.shp_beta = self.param_list['shp_beta'] = shp_beta
        self.rte_beta = self.param_list['rte_beta'] = rte_beta
        

        # State variables
        self.core_Q[:] = 1.
        self.beta = 1.
        self.gam = 1.
        self.Theta_K = np.ones(self.n_subpops)
        self.Theta_IK = np.ones((self.n_subjects, self.n_subpops))
        self.Theta_SK = np.ones((self.n_samples, self.n_subpops))
        self.H_K = np.zeros(self.n_subpops, dtype=int)
        self.H_IK = np.zeros((self.n_subjects, self.n_subpops), dtype=int)
        self.Z_LK = np.zeros((self.n_loci, self.n_subpops), dtype=np.int32)
        self.H_T = np.zeros((self.n_loci, self.n_subpops), dtype='double')
        # print(self.H_T.flags)
        self.phi_lk = np.sum(self.H_T, axis=0)

        # Auxiliary structures
        self.P_XLG = np.zeros((self.n_threads, self.n_loci, self.n_genotypes))

    cdef list _get_variables(self):
        """
        Return variable names, values, and sampling methods for testing.

        MUST BE IN TOPOLOGICAL ORDER!
        """
        variables = [('beta', self.beta, self._update_beta),
                     ('gam', self.gam, self._update_gam),
                     # ('core_Q', self.core_Q, self._update_core_Q),
                     # ('b_M', self.b_M, self._update_b_M),
                     ('Theta_K', self.Theta_K, self._update_Theta_tree_),
                     ('Theta_IK', self.Theta_IK, self._dummy_update),
                     ('Theta_SK', self.Theta_SK, self._dummy_update),
                     ('Z_LK', self.Z_LK, self._update_Z_LK),
                     # ('H_T', self.H_T, self._dummy_update),
                     # ('phi_lk', self.phi_lk, self._dummy_update),
                     ('mtx_MKD', self.mtx_MKD, self._update_mtx_MKD),
                     ('Y_MKD', self.Y_MKD, self._update_Y_PQ),
                     ('Y_Q', self.Y_Q, self._dummy_update),
                     ('H_IK', self.H_IK, self._update_H_IK_and_H_K),
                     ('H_K', self.H_K, self._dummy_update)]
        return variables
    
    def set_state(self, state):
        for key, var, _ in self._get_variables():
            if key in state.keys():
                state_var = state[key]
                if key == 'gam':
                    self.gam = state_var
                elif key == 'beta':
                    self.beta = state_var    
                else:
                    assert var.shape == state_var.shape
                    for idx in np.ndindex(var.shape):
                        var[idx] = state_var[idx]
        self._update_cache()

    cdef void _initialize_state(self, dict state={}):
        """
        Initialize internal state.
        """
        for key, var, update_func in self._get_variables():
            if key in state.keys():
                state_var = state[key]
                if key == 'gam':
                    self.gam = state_var
                elif key == 'beta':
                    self.beta = state_var    
                else:
                    assert var.shape == state_var.shape
                    for idx in np.ndindex(var.shape):
                        var[idx] = state_var[idx]
            else:
                output = update_func(self, update_mode=self._INITIALIZE_MODE)
                exit_if(output, 'updating %s' % key)

    def generate(self):
        self._generate_state()
        self._generate_data()
        data_I_LB = self.create_data_set()
        state = dict(self.get_state())
        return data_I_LB, state

    def create_data_set(self):

        data_SM = np.zeros((self.n_samples, self.n_loci * self.n_bases), dtype=long)

        for p in range(self.n_nonzero):
            y = self.nonzero_data_P[p]
            s, m = self.nonzero_subs_PM[p]
            data_SM[s, m] = y

        data_SLB = data_SM.reshape((self.n_samples, self.n_loci, self.n_bases))

        s = 0
        data_I_LB = []
        for Ji in self.n_biopsies_I:
            data_I_LB.append(data_SLB[s:s+Ji])
            s += Ji

        return data_I_LB

    def fit(self, data, n_itns=1000, initialize=True, verbose=1,
            impute_after=0, schedule={}, fix_state={}, init_state={}):

        if len(data) == self.n_samples:
            assert isinstance(data, np.ndarray)
            assert data.shape == (self.n_samples, self.n_loci * self.n_bases)
            data_SM = data
        else:
            assert len(data) == self.n_subjects
            assert all(data[i].shape == (self.n_biopsies_I[i], self.n_loci, self.n_bases)
                for i in range(self.n_subjects))

            s = 0
            data_SM = np.zeros((self.n_samples, self.n_loci * self.n_bases), dtype=long)
            for i, Ji in enumerate(self.n_biopsies_I):
                data_SM[s:s+Ji] = data[i].reshape((Ji, -1))
                s += Ji

        assert data_SM.shape == (self.n_samples, self.n_loci * self.n_bases)
        assert data_SM.dtype == long

        # print('iter')
        # for l in range(self.n_loci):
        #     for g in range(self.n_genotypes):
        #         if g == 0:
        #             self.T_GB[l, g, 0] = 0
        #         elif g==1:
        #             self.T_GB[l, g, 0] = np.mean(data_SM[:, l])
        #         else:
        #             self.T_GB[l, g, 0] = np.max(data_SM[:, l])

        # print(self.Y_MKD.shape)
        # print('0', np.asarray(self.Y_MKD[0]))
        # print('1', np.asarray(self.Y_MKD[1]))
        # print('theta_sk', np.asarray(self.Theta_SK)[:,1])
        # print(self.Y_XMKD.shape)

        super(GaP_HDP_Seq, self).fit(data=data_SM,
                                     n_itns=n_itns,
                                     initialize=initialize,
                                     verbose=verbose,
                                     impute_after=impute_after,
                                     schedule=schedule,
                                     fix_state=fix_state,
                                     init_state=init_state)

    cdef void _generate_state(self):
        """
        Generate internal state.
        """
        for key, _, update_func in self._get_variables():
            if key not in ['Y_MKD', 'Y_Q', 'H_IK', 'H_K']:
                update_func(self, update_mode=self._GENERATE_MODE)

    cdef void _generate_data(self):
        self._update_Y_PQ(update_mode=self._GENERATE_MODE)
        self._update_H_IK_and_H_K(update_mode=self._GENERATE_MODE)

    cdef void _update_core_Q(self, int update_mode):
        self.core_Q[:] = 1.

    cdef void _update_mtx_MKD(self, int update_mode):
        # self._update_mtx_m_KD(1, update_mode)
        pass

    cdef void _update_b_M(self, int update_mode):
        """There are no gamma-distributed factors in this model.

        This overwrites the method from apf.pyx that updates the
        rate hyperprior for any gamma-distributed factors.
        """
        pass

    cdef void _update_Z_LK(self, int update_mode): ###### change ##########
        cdef:
            np.npy_intp tid, k, g, l, m, b, i
            double lse, count
            gsl_rng * rng

        for k in prange(self.n_subpops, schedule='static', nogil=True):
            tid = self._get_thread()
            rng = self.rngs[tid]

            for g in range(self.n_genotypes):
                for l in range(self.n_loci):
                    self.P_XLG[tid, l, g] = log(self.A_LG[l, g])

                if update_mode == self._INFER_MODE:
                    m = -1
                    for l in range(self.n_loci):
                        for b in range(self.n_bases):
                            m = m + 1
                            count = 0.
                            
                            for i in range(self.n_samples):
                                count = count + self.Theta_SK[i, k]
                            self.P_XLG[tid, l, g] += _logpmf_poisson(self.Y_MKD[1, k, m], count * self.T_GB[l, g, b]) ##changed 2021-12-20-15-51 l-> m
                            # self.P_XLG[tid, l, g] += _logpmf_poisson(self.Y_MKD[1, k, l], count * self.T_GB[g, b])
            for l in range(self.n_loci):
            	# if self.P_XLG[tid, l] = []
                # with gil:
                #     print(self.P_XLG.shape, np.array(self.P_XLG[tid,l]), 'l', l, 'k', k)
                #     print(_sample_logcategorical(rng, self.P_XLG[tid, l]))
                self.Z_LK[l, k] = _sample_logcategorical(rng, self.P_XLG[tid, l])
                self.H_T[l, k] = self.T_GB[l, self.Z_LK[l,k], 0]
                # if self.Z_LK[l, k] == 0:
                #     self.H_T[l, k] = self.Z_LK[l, k] * 1
                # elif self.Z_LK[l, k] == 1:
                #     self.H_T[l, k] = self.Z_LK[l, k] * 100
                # else:
                #     self.H_T[l, k] = self.Z_LK[l, k] * 500     
                
                if self.Z_LK[l, k] < 0:
                    with gil:
                        print(np.array(self.P_XLG[tid, l]))

            m = -1
            for l in range(self.n_loci):
                for b in range(self.n_bases):
                    m = m + 1
                    self.mtx_MKD[1, k, m] = self.T_GB[l, self.Z_LK[l, k], b]
                    # self.mtx_MKD[1, k, m] = self.T_GB[self.Z_LK[l, k], b]
        # self.H_T = self.Z_LK
        self.phi_lk = np.sum(self.H_T, axis=0)
	
	

    cdef void _update_beta(self, int update_mode): ### gamma ####
        cdef:
            double shp, rte

        shp = rte = 1.

        if update_mode == self._INFER_MODE:
            shp += self.gam
            rte += np.sum(self.Theta_K)

        self.beta = _sample_gamma(self.rng, shp, 1./rte)

    cdef void _update_gam(self, int update_mode): #### rho_0 ####
        cdef:
            np.npy_intp i, k
            double shp, rte, tmp

        shp, rte = self.shp_gam, self.rte_gam
        # phi_lk = np.sum(self.H_T, axis=1)

        if update_mode == self._INFER_MODE:
            for k in prange(self.n_subpops, schedule='static', nogil=True):
                rng = self.rngs[self._get_thread()]
                shp += _sample_crt(rng, self.H_K[k], self.gam / self.n_subpops)

                tmp = 0
                for i in range(self.n_subjects):
                    # tmp += log1p(self.n_biopsies_I[i] * log1p(self.n_loci)) ###### change ##########
                    tmp += log1p(self.n_biopsies_I[i] * log1p(self.phi_lk[k]))
                rte += log1p(tmp / self.beta) / self.n_subpops

        self.gam = _sample_gamma(self.rng, shp, 1./rte)

    cdef void _update_Theta_tree_(self, int update_mode):
        cdef:
            np.npy_intp k, i, j, s
            double prior_rte, shp_k, rte_k, shp_ik, rte_ik, shp_sk, rte_sk
            gsl_rng * rng

        prior_rte = 1.
        # phi_lk = np.sum(self.H_T, axis=1)

        for k in prange(self.n_subpops, schedule='static', nogil=True): # k
            rng = self.rngs[self._get_thread()]

            shp_k = self.gam / float(self.n_subpops)
            rte_k = self.beta

            if update_mode == self._INFER_MODE:
                shp_k = shp_k + self.H_K[k]
                for i in range(self.n_subjects):
                    # rte_k = rte_k + log1p(self.n_biopsies_I[i] * log1p(self.n_loci)) ### change ####
                    rte_k = rte_k + log1p(self.n_biopsies_I[i] * log1p(self.phi_lk[k]))

            self.Theta_K[k] = _sample_gamma(rng, shp_k, 1./rte_k)

            s = -1
            for i in range(self.n_subjects): #### i, k #########

                shp_ik = self.Theta_K[k]
                rte_ik = prior_rte

                if update_mode == self._INFER_MODE:
                    shp_ik = shp_ik + self.H_IK[i, k]  #### H_IK[i, k] -> w_ik
                    # rte_ik = rte_ik + self.n_biopsies_I[i] * log1p(self.n_loci) ### change ####
                    rte_ik = rte_ik + self.n_biopsies_I[i] * log1p(self.phi_lk[k])

                self.Theta_IK[i, k] = _sample_gamma(rng, shp_ik, 1./rte_ik)

                for j in range(self.n_biopsies_I[i]):
                    s = s + 1

                    shp_sk = self.Theta_IK[i, k]
                    rte_sk = prior_rte

                    if update_mode == self._INFER_MODE: ###i j k
                        shp_sk = shp_sk + self.Y_MKD[0, k, s] #### yijk -> Y_MKD 
                        # rte_sk = rte_sk + self.n_loci #### change ####
                        rte_sk = rte_sk + self.phi_lk[k]

                    self.Theta_SK[s, k] = _sample_gamma(rng, shp_sk, 1./rte_sk)
                    self.mtx_MKD[0, k, s] = self.Theta_SK[s, k]

    cdef void _update_H_IK_and_H_K(self, int update_mode): #### do not need change ######### H_IK = wik', H_k = w_k''
        cdef:
            np.npy_intp k, i, j, s
            int m_sk, m_ik
            double r_sk, r_ik
            gsl_rng * rng

        self.H_K[:] = 0
        self.H_IK[:] = 0

        for k in prange(self.n_subpops, schedule='static', nogil=True):
            rng = self.rngs[self._get_thread()]

            s = -1
            for i in range(self.n_subjects):
                for j in range(self.n_biopsies_I[i]):
                    s = s + 1

                    m_sk = self.Y_MKD[0, k, s]
                    r_sk = self.Theta_IK[i, k]

                    if m_sk > 0:
                        self.H_IK[i, k] += _sample_crt(rng, m_sk, r_sk) ### m_sk = yij(k) r_sk = theta_ik'

                m_ik = self.H_IK[i, k] ### w_ik
                r_ik = self.Theta_K[k] ### theta_k''

                if m_ik > 0:
                    self.H_K[k] += _sample_crt(rng, m_ik, r_ik)

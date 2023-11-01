import scipy.stats as st
import numpy.random as rn
import numpy as np
from apf.models.gap_hdpseq import GaP_HDP_Seq
import sys
from path import Path
sys.path.append(Path(__file__).parent.parent)


def main():
    """Main method."""
    shp_gam, rte_gam = 10, 0.1
    shp_beta, rte_beta = 10, 10

    n_loci = 7
    n_biopsies_I = [5, 5, 6]
    n_subpops = 4

    n_threads = 2
    # mask_p = 0.15
    seed = rn.randint(10000)

    # seed = 1586
    rn.seed(seed)
    print('seed: %d' % seed)

    model = GaP_HDP_Seq(n_loci=n_loci,
                        n_biopsies_I=n_biopsies_I,
                        n_subpops=n_subpops,
                        shp_gam=shp_gam,
                        rte_gam=rte_gam,
                        shp_beta=shp_beta,
                        rte_beta=rte_beta,
                        seed=seed,
                        n_threads=n_threads)

    # Y = np.zeros(data_shp, dtype=np.int32)
    # mask = None
    # if mask_p > 0:
    #     mask = rn.binomial(1, mask_p, size=data_shp)
    #     percent_missing = np.ceil(100 * mask.sum() / float(mask.size))
    #     print('%d%% missing' % percent_missing)

    # data = np.ma.array(Y, mask=mask)
    # model._initialize_data(data)

    def get_schedule_func(burnin=0, update_every=1):
        return lambda x: x >= burnin and x % update_every == 0

    schedule = {'gam': get_schedule_func(0, 1),
                'beta': get_schedule_func(0, 1),
                'core_Q': get_schedule_func(np.inf, 1),
                'b_M': get_schedule_func(np.inf, 1),
                'mtx_MKD': get_schedule_func(np.inf, 1),
                'Theta_K': get_schedule_func(0, 1),
                'Theta_IK': get_schedule_func(0, 1),
                'Theta_SK': get_schedule_func(0, 1),
                'Z_LK': get_schedule_func(0, 1),
                'Y_MKD': get_schedule_func(0, 1),
                'Y_Q': get_schedule_func(0, 1),
                'H_IK': get_schedule_func(0, 1),
                'H_K': get_schedule_func(0, 1)}

    var_funcs = {}
    model.schein(6000, var_funcs=var_funcs, schedule=schedule)


if __name__ == '__main__':
    main()

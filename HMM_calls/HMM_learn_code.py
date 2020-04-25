import numpy as np
from hmmlearn import hmm
np.random.seed(42)

table_1 = np.loadtxt("~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/HMM_calling/features_geneNORMcolNrpe_antiSense_methDiff_BothContainConvert.bed")

remodel_1 = hmm.GaussianHMM(n_components=4, covariance_type="full", n_iter=100)

remodel_1.fit(table_1)


Z1 = remodel_1.predict(table_1)

np.savetxt('~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/HMM_calling/hmm_1_output_full_3feature_4st_op3_geneNORM_BothContainConvert.txt', Z1, fmt='%i')

Z1_prob = remodel_1.predict_proba(table_1)

np.savetxt('~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/HMM_calling/option3_geneNORM_BothContainConvert_probability_4feature_1.txt', Z1_prob, fmt='%1.3f')

Z1_trans = remodel_1.transmat_

np.savetxt('~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/HMM_calling/transition_matrix_values.txt', Z1_trans, fmt='%1.3f')

Z1_LLE = remodel_1._compute_log_likelihood(table_1)

np.savetxt('~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/HMM_calling/log_likelihood_enission_values.txt', Z1_LLE, fmt='%1.3f')

Z_score = remodel_1.score(table_1)
Z_score
#2212479.3037514784

Z_decode = remodel_1.decode(table_1)
Z_decode
#(2209006.6104485556, array([0, 0, 0, ..., 0, 0, 0]))

exit()


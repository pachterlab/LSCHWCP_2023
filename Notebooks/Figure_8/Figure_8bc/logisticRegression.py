#!/usr/bin/env python
# coding: utf-8

# Logistic regression of host gene expression to predict viral presence
# last modified: May 4th, 2024

import numpy as np
import pandas as pd
import pickle
import sys
from io import StringIO

import scipy
from sklearn.linear_model import LogisticRegression

import anndata
import argparse

# Turn off warnings
import warnings

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()
parser.add_argument("-covariates", "--covariates_kind", default="none", type=str)
parser.add_argument("-regularization", "--regularization", default="l1", type=str)
parser.add_argument("-genes_kind", "--genes_kind", default="hv", type=str)
parser.add_argument("-viruses_kind", "--viruses_kind", default="mac", type=str)
parser.add_argument("-scramble", "--scramble", default=False, type=bool)
parser.add_argument("-control", "--control", default="nocontrol", type=str)
parser.add_argument("-matrix", "--matrix", default="fullM", type=str)
parser.add_argument("-random_seed", "--random_seed", default=0, type=int)

args = parser.parse_args()
covariates_kind = args.covariates_kind
regularization = args.regularization
genes_kind = args.genes_kind
viruses_kind = args.viruses_kind
scramble = args.scramble
control = args.control
matrix = args.matrix
RANDOM_SEED = args.random_seed

np.random.seed(RANDOM_SEED)


def nd(array):
    return np.array(array).flatten()


# Load data
virus_data_path = "virus_no_mask.h5ad"
host_data_path = "macaque_QC_norm_leiden_celltypes.h5ad"

virus_adata = anndata.read_h5ad(virus_data_path)
host_adata = anndata.read_h5ad(host_data_path)

# Add unique barcode column
host_adata.obs["unique_bc"] = host_adata.obs["sample_barcode"].astype(
    "str"
) + host_adata.obs["barcode"].astype("str")
virus_adata.obs["unique_bc"] = virus_adata.obs["sample_barcode"].astype(
    "str"
) + virus_adata.obs.index.astype("str")
host_adata.obs.index = host_adata.obs["unique_bc"]
virus_adata.obs.index = virus_adata.obs["unique_bc"]

# Remove non macaque genes
host_adata = host_adata[
    host_adata.obs["species"] == "macaca_mulatta",
    host_adata.var["species"] == "macaca_mulatta",
]
# Remove null cell types
virus_adata = virus_adata[virus_adata.obs["celltype"].notnull(), :]

# Filter the host anndata matrix to contain only cells in filtered virus adata
host_adata = host_adata[host_adata.obs.unique_bc.isin(virus_adata.obs.unique_bc), :]

# Subset for top half of cells based on raw expression
if matrix == "halfM":
    summed_raw_counts = np.array(host_adata.raw.X.sum(axis=1)).flatten()
    thresh = np.quantile(summed_raw_counts, 0.5)
    host_adata = host_adata[summed_raw_counts > thresh]
    virus_adata = virus_adata[summed_raw_counts > thresh]

# Filter the host anndata matrix to only contain macaque genes and the viral cells
if (
    genes_kind == "all"
):  # options: 'all', 'hv', 'threshN' with N being the lowest count sum over all cells to keep a gene
    host_adata = host_adata[host_adata.obs.unique_bc.isin(virus_adata.obs.unique_bc), :]
elif genes_kind == "hv":
    host_adata = host_adata[
        host_adata.obs.unique_bc.isin(virus_adata.obs.unique_bc),
        host_adata.var.highly_variable == True,
    ]
elif "thresh" in genes_kind:
    thresh = int(genes_kind[6:])
    host_adata = host_adata[
        host_adata.obs.unique_bc.isin(virus_adata.obs.unique_bc),
        host_adata.X.sum(axis=0) >= thresh,
    ]

# Train and test using only neutrophils
if control == "neutrophils":
    index = virus_adata.obs["celltype"].isin(["Immature neutrophils"])
    virus_adata = virus_adata[index, :]
    host_adata = host_adata[index, :]
    print("USING NEUTROPHIL CONTROL")

# Select test cells such that they are of the same cell types as training cells
if control == "equalprop":
    print("USING EQUAL PROPORTIONS")

# Define viruses to train models on
# All 'macaque only' and 'shared' viruses
top_viruses_supp = [
    "u39566",
    "u102540",
    "u11150",
    "u10",
    "u288819",
    "u290519",
    "u10240",
    "u183255",
    "u1001",
    "u100291",
    "u103829",
    "u110641",
    "u181379",
    "u202260",
    "u135858",
    "u101227",
    "u100188",
    "u27694",
    "u34159",
    "u100245",
    "u10015",
    "u100733",
    "u100173",
    "u100196",
    "u100599",
    "u100644",
    "u100296",
    "u100017",
    "u100002",
    "u100012",
    "u100024",
    "u100048",
    "u100302",
    "u100074",
    "u100289",
    "u100026",
    "u100111",
    "u100139",
    "u100154",
    "u100251",
    "u100177",
    "u100215",
    "u100049",
    "u100000",
    "u100001",
    "u100007",
    "u100004",
    "u100011",
    "u100093",
    "u100116",
    "u100019",
    "u100076",
    "u100028",
    "u100153",
    "u100031",
    "u100145",
    "u102324",
    "u134800",
]
# Viruses shown in manuscript main figures
top_viruses_top7 = [
    "u10",
    "u102540",
    "u11150",
    "u202260",
    "u39566",
    "u134800",
    "u102324",
]

if viruses_kind == "supp":
    top_viruses = top_viruses_supp
    print("Using all 'macaque only' and 'shared' viruses")
elif viruses_kind == "top7":
    top_viruses = top_viruses_top7
    print("Using 7 viruses shown in manuscript main figures")

# create smaller adata containing only the viruses of interest
virus_adata = virus_adata[:, virus_adata.var.index.isin(top_viruses)]
clusters = np.array(virus_adata.obs.celltype_clusters.unique())

# order
top_viruses = virus_adata.var.index.tolist()


# function for data selection
def extract_training_data_new(
    host_adata, virus_adata_, covariate_array=None, equal_pos_neg=True, return_idx=True
):
    """Returns training data for a logistic regression model.
    This time, only cells in top 3rd will be considered for training.

    parameters
    ------
    host_adata : contains gene expression matrix
    virus_adata : contains viral counts
    virus : which virus to extract
    covariates : list of covariates to include, must be column names in the adata.obs dataframe

    returns
    ------
    X : numpy array of shape sample by genes
    Y : numpy array shape sample by 1, where 0 or 1 is presence of absence of viruses in the sample
    """

    virus_high = virus_adata_[virus_adata_.X.sum(axis=1) > 0]
    virus_low = virus_adata_[virus_adata_.X.sum(axis=1) == 0]

    # equal positive and negative virus samples
    if equal_pos_neg == True:
        num_high = int(len(virus_high) / 2)
        num_low = int(len(virus_low) / 2)

        if num_high > num_low:
            N = num_low
        elif num_low >= num_high:
            N = num_high

        if control == "equalprop":
            virus_high = virus_high[
                np.random.choice(np.arange(num_high), size=N, replace=False)
            ]

            cts = virus_high.obs["celltype"].unique()

            n_high_ct = len(virus_high[virus_high.obs["celltype"] == cts[0], :])
            virus_low_ = virus_low[virus_low.obs["celltype"] == cts[0], :]
            virus_low_ = virus_low_[
                np.random.choice(
                    np.arange(len(virus_low_)), size=n_high_ct, replace=False
                )
            ]

            for ct in cts[1:]:
                n_high_ct = len(virus_high[virus_high.obs["celltype"] == ct, :])
                virus_low_ct = virus_low[virus_low.obs["celltype"] == ct, :]
                virus_low_ = virus_low_.concatenate(
                    virus_low_ct[
                        np.random.choice(
                            np.arange(len(virus_low_ct)), size=n_high_ct, replace=False
                        )
                    ]
                )

            virus_low = virus_low_

        else:
            virus_high = virus_high[
                np.random.choice(np.arange(num_high), size=N, replace=False)
            ]
            virus_low = virus_low[
                np.random.choice(np.arange(num_low), size=N, replace=False)
            ]

    Y_high = virus_high.X
    X_high = host_adata[host_adata.obs.unique_bc.isin(virus_high.obs.unique_bc), :].X

    Y_low = virus_low.X
    X_low = host_adata[host_adata.obs.unique_bc.isin(virus_low.obs.unique_bc), :].X

    X = scipy.sparse.vstack((X_high, X_low))
    Y = scipy.sparse.vstack((Y_high, Y_low))
    virus_high_low = virus_high.concatenate((virus_low))

    # reset index, which changed after concatenation
    virus_high_low.obs.index = (
        virus_high.obs.unique_bc.values.tolist()
        + virus_low.obs.unique_bc.values.tolist()
    )

    if covariate_array != None:
        index = virus_adata.obs.unique_bc.isin(virus_high_low.obs.unique_bc)
        covariate_array_subset = covariate_array[index, :]
        #         covariate_array = one_hot_covariates(virus_high_low,covariate_list,sparse=True)
        X = scipy.sparse.hstack((X, covariate_array_subset))

    if return_idx:
        return (
            np.array(X.todense()),
            np.array(Y.todense()),
            virus_high_low.obs.unique_bc,
        )
    else:
        return (np.array(X.todense()), np.array(Y.todense()))


# one hot encode the covariates
def one_hot_covariates(adata, covariate_list, return_sparse=True):
    """One hot encodes the covariates in covariate_list."""

    covariate_array = np.array(pd.get_dummies(adata.obs[covariate_list[0]]))

    for cov in covariate_list[1:]:
        covariate_array = np.concatenate(
            (covariate_array, np.array(pd.get_dummies(adata.obs[cov]))), axis=1
        )

    if return_sparse:
        return scipy.sparse.csr_array(covariate_array)
    else:
        return covariate_array


# deal with covariates
if covariates_kind == "none":
    covariates_array = None
    num_covariates = 0

elif covariates_kind == "donor_time":
    covariates_array = one_hot_covariates(
        host_adata, ["donor_animal", "dpi_clean_merged"]
    )
    num_covariates = covariates_array.shape[1]

num_genes = host_adata.shape[1]

result_dict = {
    "num_training": np.ones(len(top_viruses)),
    "test_sensitivity": np.ones(len(top_viruses)),
    "test_specificity": np.ones(len(top_viruses)),
    "test_score": np.ones(len(top_viruses)),
    "num_training": np.ones(len(top_viruses)),
    "training_idx": [],
    "test_score": np.ones(len(top_viruses)),
    "weights": np.ones((len(top_viruses), num_genes + num_covariates + 1)),
    "true_all": [],
    "probabilities_all": [],
    "loss": [],
}


for i, v in enumerate(top_viruses):
    print(str(v) + viruses_kind)
    X, Y, training_idx = extract_training_data_new(
        host_adata, virus_adata[:, v], covariate_array=covariates_array, return_idx=True
    )

    print(len(Y), "length of X, and sum:,", np.sum(np.array(Y)))
    result_dict["training_idx"].append(training_idx)

    # randomly shuffle training data
    rand_index = np.arange(len(X))
    np.random.shuffle(rand_index)
    X = X[rand_index, :]
    Y = Y[rand_index, :]

    if scramble == True:
        np.random.shuffle(rand_index)
        Y = Y[rand_index, :]

    result_dict["num_training"][i] = len(X)

    # Record std out to catch loss
    old_stdout = sys.stdout
    sys.stdout = mystdout = StringIO()

    # train model
    print(len(Y), "length of X, and sum:,", np.sum(np.array(Y)))
    model = LogisticRegression(
        random_state=RANDOM_SEED, penalty=regularization, solver="liblinear"
    ).fit(X, Y)

    # Record loss
    sys.stdout = old_stdout
    loss_history = mystdout.getvalue()
    loss_list = []
    for line in loss_history.split("\n"):
        if len(line.split("loss: ")) == 1:
            continue
        loss_list.append(float(line.split("loss: ")[-1]))
    result_dict["loss"].append(loss_list)

    # test model
    X_test = host_adata[~host_adata.obs.unique_bc.isin(training_idx), :].X
    if covariates_kind != "none":
        covariates_test = covariates_array[
            ~host_adata.obs.unique_bc.isin(training_idx), :
        ]
        X_test = np.array(scipy.sparse.hstack((X_test, covariates_test)).todense())
    Y_test = np.array(
        virus_adata[~host_adata.obs.unique_bc.isin(training_idx), v].X.todense()
    ).flatten()

    # score
    result_dict["test_score"][i] = model.score(X_test, Y_test)
    result_dict["probabilities_all"].append(model.predict_proba(X_test).flatten())
    result_dict["true_all"].append(np.array(Y_test).flatten())

    # sensitivity -- true positives
    virus_adata_test = virus_adata[~virus_adata.obs.unique_bc.isin(training_idx), v]
    virus_adata_test = virus_adata_test[virus_adata_test.X.sum(axis=1) > 0]

    X_test = np.array(
        host_adata[
            host_adata.obs.unique_bc.isin(virus_adata_test.obs.unique_bc), :
        ].X.todense()
    )
    if covariates_kind != "none":
        covariates_test = covariates_array[
            host_adata.obs.unique_bc.isin(virus_adata_test.obs.unique_bc), :
        ]
        X_test = np.array(scipy.sparse.hstack((X_test, covariates_test)).todense())
    Y_test = np.array(virus_adata_test.X.todense()).flatten()

    result_dict["test_sensitivity"][i] = model.score(X_test, Y_test)

    # specificity -- true negatives
    virus_adata_test = virus_adata[~virus_adata.obs.unique_bc.isin(training_idx), v]
    virus_adata_test = virus_adata_test[virus_adata_test.X.sum(axis=1) == 0]

    X_test = np.array(
        host_adata[
            host_adata.obs.unique_bc.isin(virus_adata_test.obs.unique_bc), :
        ].X.todense()
    )
    if covariates_kind != "none":
        covariates_test = covariates_array[
            host_adata.obs.unique_bc.isin(virus_adata_test.obs.unique_bc), :
        ]
        X_test = np.array(scipy.sparse.hstack((X_test, covariates_test)).todense())
    Y_test = np.array(virus_adata_test.X.todense())

    result_dict["test_specificity"][i] = model.score(X_test, Y_test)

    result_dict["weights"][i, :] = np.hstack(
        (model.intercept_[:, None], model.coef_)
    ).flatten()

# save results
result_dict["viruses"] = top_viruses

if scramble == True:
    with open(
        f"{viruses_kind}_viruses_{genes_kind}_genes_{matrix}_cov_{covariates_kind}_{regularization}_{control}_scramble_{RANDOM_SEED}.pickle",
        "wb",
    ) as handle:
        pickle.dump(result_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
else:
    with open(
        f"{viruses_kind}_viruses_{genes_kind}_genes_{matrix}_cov_{covariates_kind}_{regularization}_{control}_{RANDOM_SEED}.pickle",
        "wb",
    ) as handle:
        pickle.dump(result_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

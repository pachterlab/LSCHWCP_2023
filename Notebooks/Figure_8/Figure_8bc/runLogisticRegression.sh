


# hv genes
python3 logisticRegression.py --covariates_kind 'donor_time' --genes_kind 'hv' --regularization 'l2' --viruses_kind "supp"

python3 logisticRegression.py --covariates_kind 'none' --genes_kind 'hv' --regularization 'l2' --viruses_kind "supp"


# hv genes scrambles
python3 logisticRegression.py --covariates_kind 'donor_time' --genes_kind 'hv' --regularization 'l2' --viruses_kind "supp" --scramble True

python3 logisticRegression.py --covariates_kind 'none' --genes_kind 'hv' --regularization 'l2' --viruses_kind "supp" --scramble True


# hv genes controls
python3 logisticRegression.py --covariates_kind 'donor_time' --genes_kind 'hv' --regularization 'l2' --viruses_kind "top7" --control "equalprop"

python3 logisticRegression.py --covariates_kind 'donor_time' --genes_kind 'hv' --regularization 'l2' --viruses_kind "top7" --control "equalprop" --scramble True




# ALL genes
python3 logisticRegression.py --covariates_kind 'donor_time' --genes_kind 'all' --regularization 'l2' --viruses_kind "supp"

python3 logisticRegression.py --covariates_kind 'none' --genes_kind 'all' --regularization 'l2' --viruses_kind "supp"


# ALL genes scramble
python3 logisticRegression.py --covariates_kind 'donor_time' --genes_kind 'all' --regularization 'l2' --viruses_kind "supp" --scramble True

python3 logisticRegression.py --covariates_kind 'none' --genes_kind 'all' --regularization 'l2' --viruses_kind "supp" --scramble True


#  ALL genes controls
python3 logisticRegression.py --covariates_kind 'donor_time' --genes_kind 'all' --regularization 'l2' --viruses_kind "top7" --control "equalprop"

python3 logisticRegression.py --covariates_kind 'donor_time' --genes_kind 'all' --regularization 'l2' --viruses_kind "top7" --control "equalprop" --scramble True




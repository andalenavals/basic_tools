## EXAMPLES

This is is set of particular examples where some of the base base
codes defiend in the basic tools ares used.

#filter_exp.py

This code filter a list of exposures, based on a
physical criterium, using information from a extra file.
For instance,

python filter_exp.py --explist 'ally3.riz' --file=mean_ellip_byexp.fits >> newout.riz

will save in newout.riz a list of exposure with mean_os_e lower than 0.11
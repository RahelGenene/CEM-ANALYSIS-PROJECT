CEM Project - HIV/ART risk prediction and clinical intelligence dashboard
Rahel Genene (1601570) and Hanna Getachw (1601241), DBU Data Science internship, MoH Ethiopia

Deployed dashboard: https://rahelgenene.shinyapps.io/cem_deploy/

dataset/         raw CEM.csv (AHRI REDCap export)
notebook/        notebooks only, run in this order
  00_data_understanding.ipynb
  01_pre_processing.ipynb
  02_EDA.ipynb
  03_ML_modeling_mortality.ipynb
  04_ML_modeling_unfavorable_outcome.ipynb
output/          images, tables, data and models, one folder per notebook
  00_data_understanding/ 01_pre_processing/ 02_EDA/ 03_ML_modeling_mortality/ 04_ML_modeling_unfavorable_outcome/
dashboard/       app.R, data file, README with the deployed link, link shortcut
documentation/   internship report (.docx)
requirements.txt  Python packages for the notebooks

Notes
- Notebooks 00 and 04 use relative paths and were run from the notebook/ folder; their outputs are saved in output/.
- Notebooks 01, 02 and 03 are the original working notebooks and still read from your local Lenovo path or the working folder. To rerun them, point them to dataset/CEM.csv and to the files in output/.
- Notebook 04 rebuilds the unfavorable-outcome model from the pre-processing outputs; its metrics differ slightly from the figures in the report (which came from the original run).
- output/01_pre_processing/figures/cd4_by_visit_chart.png is an older version, superseded by cd4_by_visit_chart_clean.png.

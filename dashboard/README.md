CEM Clinical Intelligence Dashboard (R Shiny)

Live (deployed on shinyapps.io): https://rahelgenene.shinyapps.io/cem_deploy/

Tabs: Executive Summary, Cohort Explorer, Risk Assessment, Facility Performance, Data Quality.

Files
- app.R                        the Shiny app
- CEM_final_preprocessed.csv   default data file (loaded from the working directory; a CSV can also be uploaded in the sidebar)
- Deployed_Dashboard.url       shortcut to the live dashboard

Run locally
1. install.packages(c("shiny","shinydashboard","DT","data.table","dplyr","tidyr","stringr","plotly","lubridate","rpart","randomForest","e1071","gbm","pROC"))
2. shiny::runApp("dashboard")   (from the project root)

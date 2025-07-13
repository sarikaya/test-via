set -e

cd paramCheckRun && rm -rf initialrun && cp -R ./../initialrun . && nextflow run -c ./../nextflow.config paramCheck.nf
if [ $INFORMATICS_USERNAME == "your_username" ]; then
    echo " -------------------------------------------------"
    printf "\nERROR: API credentials are missing.\n"
    echo "Please define INFORMATICS_USERNAME and INFORMATICS_APIKEY in your environment variables. Profile->Run Environments Tab->Click Edit Moderna Run Environment-> Enter variables into Environment Variables section."
    echo "If you need support for API credentials, please contact Nathaniel Reynolds (nathaniel.reynolds@modernatx.com)"
    exit 1
fi
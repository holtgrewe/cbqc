#!/bin/bash

# get parent dir of bash file
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# ---------------------------------------------------------------------------
# Interpret Input Parameters
# ---------------------------------------------------------------------------

# check that the IN_DIR and OUT_PREFIX are set
if [[ -z "${IN_DIR}" ]]; then
    echo "variable IN_DIR not set" 1>&2
    exit 1
fi
if [[ -z "${OUT_PREFIX}" ]]; then
    echo "variable OUT_PREFIX not set" 1>&2
    exit 1
fi

ID=$(basename ${IN_DIR})
echo "ID=${ID}"

# ---------------------------------------------------------------------------
# BASH setup
# ---------------------------------------------------------------------------

set -x  # show what we do
set -e  # exit after command failed

# create temporary directory, setup error handler to remove afterwards
TMPDIR=$(mktemp -d)
cleanup() { rm -rf ${TMPDIR}; }
trap cleanup EXIT INT TERM

# ---------------------------------------------------------------------------
# Generate HTML QC Report
# ---------------------------------------------------------------------------

fastqc_zips=${IN_DIR}/qc/*.zip
fastqc_args=
for x in ${fastqc_zips}; do
    fastqc_args=" ${fastqc_args} --fastqc-zip ${x}"
done

#virtualenv ${TMPDIR}/venv
#. ${TMPDIR}/venv/bin/activate
#pip install -r ${DIR}/requirements.txt

PYTHONPATH=${DIR} \
python -m cbpipeline.qc_report \
           --sta-file ${IN_DIR}/cov/${ID}_exom.sta \
           --out-html ${TMPDIR}/qcreport.html \
           ${fastqc_args} \
           --vcf-file ${IN_DIR}/vcf/${ID}.ra.rc.gt.CCDS.vcf

# ---------------------------------------------------------------------------
# Convert QC Report to PDF
# ---------------------------------------------------------------------------

cp ${TMPDIR}/qcreport.html ${OUT_PREFIX}.html

wkhtmltopdf ${TMPDIR}/qcreport.html ${TMPDIR}/qcreport.pdf

# there is a problem with wkhtmltopdf on the HTML with inline images, thus extract

mkdir -p ${TMPDIR}/fastqc && pushd ${TMPDIR}/fastqc
for f in ${IN_DIR}/qc/${ID}-*.zip ; do
    unzip ${f} 
done
for f in *; do
    wkhtmltopdf ${f}/fastqc_report.html ${f}.pdf
done
popd
wkhtmltopdf ${IN_DIR}/qm/qualimapReport.html ${TMPDIR}/qualimapReport.pdf

pdfjoin -o ${OUT_PREFIX}.pdf ${TMPDIR}/qcreport.pdf ${TMPDIR}/fastqc/*.pdf ${TMPDIR}/qualimapReport.pdf

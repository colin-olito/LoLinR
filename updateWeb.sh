#!/bin/sh
set -e

WEB_DIR=inst/web
VIG_DIR=${WEB_DIR}/vignettes

if [ ! -d "$VIG_DIR" ]
then
	mkdir -p ${VIG_DIR}
fi

Rscript -e 'devtools::build_vignettes()'

DOC_DIR=inst/doc
mv ${DOC_DIR}/*.html ${VIG_DIR}

VERSION=$(git rev-parse --short HEAD)
REMOTE_URL=$(git config --get remote.origin.url)

rm -rf ${WEB_DIR}/.git
git init ${WEB_DIR}
git -C ${WEB_DIR} checkout --orphan gh-pages
git -C ${WEB_DIR} add --all
git -C ${WEB_DIR} commit --no-verify -m "Update docs for version ${VERSION}"
git -C ${WEB_DIR} remote add origin -m "gh-pages" ${REMOTE_URL}
git -C ${WEB_DIR} push --force -u origin gh-pages

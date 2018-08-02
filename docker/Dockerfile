FROM genepattern/docker-r-2-15:0.1


#
# copy the Rprofile.site file, to set the default CRAN mirror
#   needed by installPackages.R
#
COPY R2.15_Rprofile.gp.site /usr/local/lib/R/etc/Rprofile.site

# add this line to the Rprofile.site, needed to install packages from https sites
#   options(download.file.method = "wget")
# 
RUN echo "\\noptions(download.file.method = \"wget\")\\n" >> /usr/local/lib/R/etc/Rprofile.site

#
# install R packages for the ExpressionFileCreator module 
#
COPY r.package.info /build/source/r.package.info
RUN Rscript --no-save --quiet --slave --no-restore \
    /build/source/installPackages.R \
    /build/source/r.package.info \
    /build/source/install.packages.log

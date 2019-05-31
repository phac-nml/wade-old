FROM r-base:latest

MAINTAINER Adrian Zetner "adrian.zetner@canada.ca"

RUN apt-get update && apt-get install -y \
    sudo \
    gdebi-core \
	net-tools \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    xtail \
    wget

# Download and install shiny server
RUN wget --no-verbose https://download3.rstudio.org/ubuntu-14.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb && \
    . /etc/environment && \
    chown shiny:shiny /var/lib/shiny-server

# Install some R packages
RUN R -e "options(repos = c(CRAN = 'http://cran.rstudio.com')); install.packages(c('remotes'));"
RUN R -e "source('https://bioconductor.org/biocLite.R'); biocLite('Rsamtools')"

# Install Pavian from github
# 	From 'galaxy' branch
# 	Use ref= to specify the branch to download. Removing ref= defaults to master
RUN R -e "remotes::install_github('phac-nml/wade')"

# Install GalaxyConnector from github
# 	Made some modifications to r-galaxy-connector specifically for Pavian hence
# 	the connector is being installed from the 'pavian' branch
#
# 	This will need to change to phac-nml/galaxy-connector
RUN R -e "remotes::install_github('justinband/r-galaxy-connector', ref='pavian')"

COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf
RUN rm -f /srv/shiny-server/*

# Move pavians app.R file into /srv/shiny-server/ to be launched
RUN R -e "sapply(list.files(system.file('shinyapp',package='wade'),full.names=TRUE),file.copy,to='/srv/shiny-server/',recursive=TRUE)"

COPY shiny-server.sh /usr/bin/shiny-server.sh
COPY monitor_traffic.sh /monitor_traffic.sh

EXPOSE 80

# Run it!
CMD ["/usr/bin/shiny-server.sh"]


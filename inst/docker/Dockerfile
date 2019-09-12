FROM r-base:3.5.2

MAINTAINER Adrian Zetner "adrian.zetner@canada.ca"

RUN apt-get update && apt-get install -y \
    sudo \
    gdebi-core \
	ncbi-blast+ \
	net-tools \
    pandoc \
    pandoc-citeproc \
	procps \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
	libxml2-dev \
    libssl-dev \
	xtail \
    nginx \
    wget \
    curl


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

# Install GalaxyConnector from github
# 	This will need to change to phac-nml/galaxy-connector
RUN R -e "remotes::install_github('phac-nml/r-galaxy-connector@*release')"

# Install dependencies before installing wade
# 	Upside: Removes the install time of wade when any changes occur.
# 	Downside: Must manage package installation from here..
RUN R -e "install.packages(c('devtools', 'dplyr', 'here', 'purrr', 'shiny', 'shinydashboard', 'shinyFiles', 'shinyWidgets', 'stringr'))"
RUN R -e "install.packages('rappdirs')"

RUN tag=$(curl --silent "https://api.github.com/repos/phac-nml/wade-data/releases/latest" | grep '"tag_name":' | sed -E 's/.*"([^"]+)".*/\1/'); curl -L "https://github.com/phac-nml/wade-data/archive/$tag.zip" --output $tag.zip; unzip $tag.zip; mv wade-data* /srv/shiny-server/data/

#COPY data/ /srv/shiny-server/data/
RUN chmod -R 757 /srv/shiny-server/data # There must be an easier way..

# Install wade from github
# 	From 'galaxy' branch
# 	Use ref= to specify the branch to download. Removing ref= defaults to master
RUN R -e "remotes::install_github('phac-nml/wade', ref='master')"

COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf

# Move WADE app.R file into /srv/shiny-server/ to be launched
RUN R -e "sapply(list.files(system.file('shinyapp',package='wade'),full.names=TRUE),file.copy,to='/srv/shiny-server/',recursive=TRUE)"

# Copy some files over
# 	shiny-server specifies launch instructions
# 		copy environment variables and such
#
# 	monitor_traffic kills the container after it hasn't been used in a while
ADD shiny-server.sh /usr/bin/shiny-server.sh
COPY monitor_traffic.sh /monitor_traffic.sh
ADD nginx-default /etc/nginx/sites-available/default

# Open port 80
EXPOSE 80 3838

# Run it!
CMD ["/usr/bin/shiny-server.sh"]


FROM r-base

ENV RENV_PATHS_CACHE=/renv/cache
ENV RENV_VERSION=0.15.4

RUN apt-get update; apt-get install libcurl4-openssl-dev libssl-dev
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

COPY ./renv.lock /renv/tmp/renv.lock

WORKDIR /renv/tmp

RUN R -e "renv::restore()"

RUN mkdir /asrt-beamsampling

WORKDIR /asrt-beamsampling
RUN mkdir -p Figures/figlist

ENTRYPOINT [ "R" ]
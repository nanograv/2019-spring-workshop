FROM nanograv/workshop-2019:b1f5f8cb04b3

COPY . ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}



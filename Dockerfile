FROM nanograv/workshop-2019:dc737ee8f3cd

COPY . ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}



FROM nanograv/workshop-2019@sha256:b1f5f8cb04b362b245611da3b861f9a006a9b670f163e646c87dffe45bd60f07

COPY . ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}



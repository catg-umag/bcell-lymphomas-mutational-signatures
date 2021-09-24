FROM r-base:4.1.1

ARG VERSION
ARG DEBIAN_FRONTEND=noninteractive

LABEL author="Diego Alvarez (dialvarezs@gmail.com)"
LABEL description="Contains Jupyter with IRkernel an some R packages (nnls,tidyverse,deconstructSigs mainly)"
LABEL version=${VERSION}


RUN apt-get update && apt-get upgrade -y

# install packages. use precompiled R packages to reduce container build time
RUN apt-get install -y \
		procps python3-pip \
		r-cran-cluster r-cran-cowplot r-cran-r.utils r-cran-factoextra r-cran-irkernel r-cran-nnls r-cran-tidyverse \
		r-cran-biocmanager r-bioc-bsgenome r-bioc-genomeinfodb r-bioc-variantannotation && \
	apt-get clean

# install jupyter
RUN pip3 install --no-cache jupyter jupyterlab nbconvert

COPY install_deps.r /root/install_deps.r
RUN Rscript /root/install_deps.r
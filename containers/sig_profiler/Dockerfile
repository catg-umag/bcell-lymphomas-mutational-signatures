FROM pytorch/pytorch:1.12.1-cuda11.3-cudnn8-runtime

ARG VERSION
ARG DEBIAN_FRONTEND=noninteractive

LABEL author="Diego Alvarez (dialvarezs@gmail.com)"
LABEL description="Contains SigProfilerExtractor and its dependencies"
LABEL version=${VERSION}
LABEL org.opencontainers.image.source https://github.com/catg-umag/bcell-lymphomas-mutational-signatures


RUN apt-get update && apt-get upgrade -y && apt-get clean

RUN /opt/conda/bin/pip install --no-cache-dir SigProfilerExtractor==${VERSION}

# test
RUN /opt/conda/bin/python3 -c "import SigProfilerExtractor"
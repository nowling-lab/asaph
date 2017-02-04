FROM debian:jessie

RUN apt-get update && \
    apt-get -y install \
    	    g++ \
	    gfortran \
	    libfreetype6 \
    	    libfreetype6-dev \
	    libpng-dev \
    	    pkg-config \
    	    python-dev \
	    python-numpy \
    	    python-pip \
	    python-scipy \
	    && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

COPY . /opt/asaph
WORKDIR /opt/asaph
RUN pip install --requirement /opt/asaph/requirements.txt
ENV PATH="/opt/asaph/bin:${PATH}"

ENTRYPOINT /bin/bash

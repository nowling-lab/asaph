FROM debian:jessie

RUN apt-get update && \
    apt-get -y install \
    	    g++ \
	        gfortran \
	        git \
	        libfreetype6 \
    	    libfreetype6-dev \
	        libpng-dev \
    	    pkg-config \
    	    python-dev \
	        python-numpy \
    	    python-pip \
	        python-scipy \
            python-seaborn \
	    && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /opt
RUN git clone https://github.com/sstephenson/bats.git
WORKDIR /opt/bats
RUN ./install.sh /usr/local

COPY . /opt/asaph
WORKDIR /opt/asaph
RUN pip install --requirement /opt/asaph/requirements.txt
ENV PATH="/opt/asaph/bin:${PATH}"

CMD /bin/bash

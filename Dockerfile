FROM ubuntu:16.04
RUN apt-get update
RUN apt-get install -qy --no-install-recommends \
curl \
g++ \
  make \
  gcc \
  pkg-config \
  ca-certificates \
  locales \
  libcurl4-openssl-dev \
  libssl-dev \
  libssh2-1-dev \
  libxml2-dev

RUN mkdir /tmp/downloads

RUN curl -sSL -o tmp.tar.gz --retry 10  https://www.python.org/ftp/python/2.7.6/Python-2.7.6.tgz && \
mkdir /tmp/downloads/python && \
tar -C /tmp/downloads/python/ --strip-components 1 -zxf tmp.tar.gz && \
cd /tmp/downloads/python/ && ./configure && make && make install && \
cd / && rm -f /tmp/downloads/tmp.tar.gz

RUN mkdir /tmp/downloads/pip && \
curl -sSL --retry 10 https://bootstrap.pypa.io/get-pip.py -o /tmp/downloads/pip/get-pip.py && \
cd /tmp/downloads/pip && \
python get-pip.py pip==19.1.1

RUN pip install numpy==1.8.2 pandas==0.14.1 argparse==1.2.1 path.py==5.1

COPY *.py /usr/local/bin/

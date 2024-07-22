FROM ubuntu:24.04

RUN apt-get update

RUN apt-get -y install \
  python3-pip \
  python3-numpy \
  python3-pandas \
  python3-scipy \
  python3-skbio \
  git
RUN apt-get install -y python3-sklearn

RUN git clone https://github.com/jwanglab/minnow
RUN cd minnow && mkdir incl && cd incl && git clone https://github.com/lh3/minimap2 && git clone https://github.com/attractivechaos/klib
WORKDIR "/minnow/incl/minimap2"
RUN make
WORKDIR "/minnow"
RUN make
WORKDIR "/"
RUN git clone https://github.com/jwanglab/fastat && cd fastat && mkdir -p incl && cd incl && git clone https://github.com/attractivechaos/klib && cd .. && make

ENV PATH="$PATH:/minnow:/minnow/incl/minimap2:/fastat"

ADD *.sh /
RUN mkdir -p /src
ADD src/*.py /src/

ENTRYPOINT ["bash", "classify.sh"]

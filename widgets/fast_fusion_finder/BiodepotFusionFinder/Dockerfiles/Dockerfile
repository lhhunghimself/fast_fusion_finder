ARG ALPINE_VERSION=3.15

FROM alpine:${ALPINE_VERSION} AS alpine-build-samtools

ARG SAMTOOLS_VERSION=1.13

RUN apk update && apk add build-base wget zlib-dev tar bzip2-dev xz-dev \
    curl-dev ncurses-dev ncurses-static curl-static zlib-static bzip2-static \
    nghttp2-static openssl-libs-static brotli-static

WORKDIR /build
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2
RUN tar -xf samtools-${SAMTOOLS_VERSION}.tar.bz2
WORKDIR /build/samtools-${SAMTOOLS_VERSION}
RUN ./configure && \
    make -j 4
RUN strip samtools

FROM alpine:${ALPINE_VERSION} AS alpine-build-minimap2
ARG MINIMAP2_VERSION=2.24
RUN apk add --no-cache build-base zlib-dev wget tar
WORKDIR /build
RUN wget -O minimap2.tar.bz2 https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}.tar.bz2
RUN tar -xjf minimap2.tar.bz2
WORKDIR /build/minimap2-${MINIMAP2_VERSION}
RUN make -j 4 && strip minimap2

FROM alpine:${ALPINE_VERSION} AS alpine-build-racon
ARG RACON_VERSION=1.5.0
RUN apk add --no-cache build-base zlib-dev cmake bash git linux-headers
WORKDIR /build
RUN git clone --recursive --branch ${RACON_VERSION} --depth 1 https://github.com/lbcb-sci/racon.git
WORKDIR /build/racon/build
RUN cmake -DCMAKE_BUILD_TYPE=Release .. && make -j 4 && strip bin/racon

FROM alpine:${ALPINE_VERSION}
ARG SAMTOOLS_VERSION=1.13
ARG MINIMAP2_VERSION=2.24
RUN apk add --no-cache bash perl libbz2 zlib libcurl xz-libs libstdc++ libgcc \
    && rm -rf /var/cache/apk/*
WORKDIR /
COPY --from=alpine-build-samtools /build/samtools-${SAMTOOLS_VERSION}/samtools /usr/local/bin/samtools
COPY --from=alpine-build-minimap2 /build/minimap2-${MINIMAP2_VERSION}/minimap2 /usr/local/bin/minimap2
COPY --from=alpine-build-racon    /build/racon/build/bin/racon                  /usr/local/bin/racon
ADD filterSam.pl /usr/local/bin/filterSam.pl
ADD filterEnds.pl /usr/local/bin/filterEnds.pl
ADD breakpointPileup.pl /usr/local/bin/breakpointPileup.pl
ADD processBreakpointDir.sh /usr/local/bin/processBreakpointDir.sh
ADD polishFusionConsensus.sh /usr/local/bin/polishFusionConsensus.sh
ADD runBFF.sh /usr/local/bin/runBFF.sh

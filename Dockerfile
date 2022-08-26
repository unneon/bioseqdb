FROM postgres:alpine

RUN apk update \
    && apk add autoconf automake bzip2-dev cmake g++ git make openssl postgresql-dev xz-dev zlib-dev \
    && git clone --depth 1 --single-branch --branch 1.5 --recursive "https://github.com/samtools/htslib.git" \
    && git clone --depth 1 --single-branch "https://github.com/lh3/bwa.git" \
    && cd /htslib && autoreconf -i && ./configure --host x86_64-pc-linux-musl && make CFLAGS="-g -O2 -fPIC" lib-static \
    && cd /bwa && make CFLAGS="-g -O2 -fPIC" libbwa.a
COPY . /bioseqdb
RUN cd /bioseqdb && mkdir build && cd build && cmake .. -DCMAKE_LIBRARY_PATH="/bwa;/htslib" -DCMAKE_CXX_FLAGS="-I /" && make install \
    && (echo "CREATE EXTENSION bioseqdb;" > /docker-entrypoint-initdb.d/bioseqdb.sql)

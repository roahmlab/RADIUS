# Copyright 2019 - 2022 The MathWorks, Inc.

# To specify which MATLAB release to install in the container, edit the value of the MATLAB_RELEASE argument.
# Use lower case to specify the release, for example: ARG MATLAB_RELEASE=r2021b
ARG MATLAB_RELEASE=r2022a

# When you start the build stage, this Dockerfile by default uses the Ubuntu-based matlab-deps image.
# To check the available matlab-deps images, see: https://hub.docker.com/r/mathworks/matlab-deps
FROM mathworks/matlab-deps:${MATLAB_RELEASE}

# Declare the global argument to use at the current build stage
ARG MATLAB_RELEASE

# Install mpm dependencies
RUN export DEBIAN_FRONTEND=noninteractive && apt-get update && \
    apt-get install --no-install-recommends --yes \
        wget \
        unzip \
        ca-certificates && \
    apt-get clean && apt-get autoremove

# Run mpm to install MATLAB in the target location and delete the mpm installation afterwards.
# If mpm fails to install successfully then output the logfile to the terminal, otherwise cleanup.
RUN wget -q https://www.mathworks.com/mpm/glnxa64/mpm && \ 
    chmod +x mpm && \
    ./mpm install \
        --release=${MATLAB_RELEASE} \
        --destination=/opt/matlab \
        --products MATLAB ROS_Toolbox Mapping_Toolbox Parallel_Computing_Toolbox Optimization_Toolbox Phased_Array_System_Toolbox DSP_System_Toolbox Signal_Processing_Toolbox Statistics_and_Machine_Learning_Toolbox || \
    (echo "MPM Installation Failure. See below for more information:" && cat /tmp/mathworks_root.log && false) && \
    rm -f mpm /tmp/mathworks_root.log && \
    ln -s /opt/matlab/bin/matlab /usr/local/bin/matlab

# Add "matlab" user and grant sudo permission.
RUN adduser --shell /bin/bash --disabled-password --gecos "" matlab && \
    echo "matlab ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers.d/matlab && \
    chmod 0440 /etc/sudoers.d/matlab

# One of the following 2 ways of configuring the license server to use must be
# uncommented.

ARG LICENSE_SERVER
# Specify the host and port of the machine that serves the network licenses 
# if you want to bind in the license info as an environment variable. This 
# is the preferred option for licensing. It is either possible to build with 
# something like --build-arg LICENSE_SERVER=27000@MyServerName, alternatively
# you could specify the license server directly using
#       ENV MLM_LICENSE_FILE=27000@flexlm-server-name
ENV MLM_LICENSE_FILE=$LICENSE_SERVER

# Alternatively you can put a license file into the container.
# You should fill this file out with the details of the license 
# server you want to use and uncomment the following line.
# COPY network.lic /opt/matlab/licenses/

# The following environment variables allow MathWorks to understand how this MathWorks 
# product (MATLAB Dockerfile) is being used. This information helps us make MATLAB even better. 
# Your content, and information about the content within your files, is not shared with MathWorks. 
# To opt out of this service, delete the environment variables defined in the following line. 
# See the Help Make MATLAB Even Better section in the accompanying README to learn more: 
# https://github.com/mathworks-ref-arch/matlab-dockerfile#help-make-matlab-even-better
ENV MW_DDUX_FORCE_ENABLE=true MW_CONTEXT_TAGS=MATLAB:DOCKERFILE:V1

# Set user and work directory
#USER matlab
#WORKDIR /home/matlab
#ENTRYPOINT ["matlab"]
#CMD [""]

# Fix matlab warning
RUN apt-get install -yq libcanberra-gtk-module libcanberra-gtk3-module




#
# ROS Installation
#

# Upgrade beforehand
RUN apt-get update -yq
RUN apt-get upgrade -yq

# Basic utilities
RUN apt-get install -yq git 
RUN apt-get install -yq wget 
RUN apt-get install -yq vim


# Install ROS
RUN apt-get install -yq lsb-release
RUN apt-get install -yq gnupg
RUN apt-get install -yq curl
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends tzdata
RUN sh -c 'echo "deb http://packages.ros.org/ros/ubuntu $(lsb_release -sc) main" > /etc/apt/sources.list.d/ros-latest.list'
RUN curl -s https://raw.githubusercontent.com/ros/rosdistro/master/ros.asc | apt-key add -
RUN apt-get update -yq
RUN apt-get install -yq ros-noetic-ros-base
RUN apt-get install -yq python3-catkin-tools
RUN apt-get install -yq ros-noetic-joy
RUN apt-get install -yq ros-noetic-tf*
RUN apt-get install -yq ros-noetic-tf2*
RUN apt-get install -yq ros-noetic-jsk-recognition-msgs
RUN apt-get install -yq libxmlrpcpp-dev librosconsole-dev
RUN apt-get install -yq g++-9
RUN apt-get install -yq clang++-12
RUN apt-get install -yq clang-format-12
RUN apt-get install -yq doxygen
RUN apt-get install -yq ripgrep

# Install newer cmake version
RUN wget https://github.com/Kitware/CMake/releases/download/v3.22.5/cmake-3.22.5-linux-x86_64.sh \
      -q -O /tmp/cmake-install.sh \
    && chmod u+x /tmp/cmake-install.sh \
    && mkdir -p /usr/bin/cmake_v3.22.5 \
    && /tmp/cmake-install.sh --skip-license --prefix=/usr/bin/cmake_v3.22.5 \
    && rm /tmp/cmake-install.sh
ENV PATH="/usr/bin/cmake/bin:${PATH}"

RUN sh -c 'echo "deb http://packages.ros.org/ros2/ubuntu $(lsb_release -sc) main" > /etc/apt/sources.list.d/ros2-latest.list'
RUN apt-get update -yq
RUN apt-get install -yq --no-install-recommends ros-foxy-ros-core
RUN apt-get install -yq ros-foxy-demo-nodes-cpp-native
RUN apt-get install -yq --no-install-recommends ros-foxy-desktop

ENV ROS1_DISTRO noetic
ENV ROS2_DISTRO foxy

RUN apt-get install -yq software-properties-common
RUN apt-get update
RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt-get -yq install python3.9 python3.9-venv
RUN apt-get install -yq python3-empy












#
# Install IPOPT
#
COPY util/coinhsl-2021.05.05.tar.gz /tmp/coinhsl-2021.05.05.tar.gz
#COPY setup-ipopt-and-hsl.sh /tmp/setup-ipopt-and-hsl.sh
#RUN /tmp/setup-ipopt-and-hsl.sh
#RUN chmod -R 777 /installs

# Ipopt Version to install, note that the specified version
# must also have a corresponding release tag in Ipopt's git repo
ARG IPOPT_VERS=3.14.4

# Directory to make new folder in
ARG TOP_LEVEL_SRC_DIR=/installs

# Location of CoinHSL tar.gz file
ARG COINHSL_TAR_GZ_DIR=/tmp

# Name of CoinHSL file without the extension
ARG COINHSL_TAR_GZ_FBASE_NAME=coinhsl-2021.05.05

# -----------------------------------------------------------------------------
# The following directories will be created
#   (parent) TOP_LEVEL_SRC_DIR/Ipopt
#   (hsl)    TOP_LEVEL_SRC_DIR/Ipopt/ThirdParty-HSL
#   (ipopt)  TOP_LEVEL_SRC_DIR/Ipopt/Ipopt-$IPOPT_VERS
# 
# HSL will be copied from COINHSL_TAR_GZ_FPATH/COINHSL_TAR_GZ_FBASE_NAME.tar.gz
#   and extracted to
#   to TOP_LEVEL_SRC_DIR/Ipopt/ThirdParty-HSL/coinhsl
# -----------------------------------------------------------------------------

# Setup directories and filepaths
ARG IPOPT_PARENT_DIR=$TOP_LEVEL_SRC_DIR/Ipopt
ARG IPOPT_WITH_VERS=Ipopt-$IPOPT_VERS
ARG COINHSL_TAR_GZ_FNAME=$COINHSL_TAR_GZ_FBASE_NAME.tar.gz
ARG COINHSL_TAR_GZ_FPATH=$COINHSL_TAR_GZ_DIR/$COINHSL_TAR_GZ_FNAME
RUN apt-get install -yq \
	g++ \
	gcc \
	gfortran \
	git \
	libatlas-base-dev \
	liblapack-dev \
	libmetis-dev \
	patch \
	pkg-config \
	unzip \
	wget

RUN mkdir -p $IPOPT_PARENT_DIR
WORKDIR $IPOPT_PARENT_DIR

# -----------
# Install HSL
# -----------

# Clone the repo
RUN git clone https://www.github.com/coin-or-tools/ThirdParty-HSL.git
WORKDIR ThirdParty-HSL

# Extract HSL to local directory
RUN tar -xvzf $COINHSL_TAR_GZ_FPATH

# Rename extracted directory for installation purposes
RUN mv $COINHSL_TAR_GZ_FBASE_NAME coinhsl

RUN ./configure && \
  make -j && \
  make install

# -------------
# Install IPOPT
# -------------

# Clone IPOPT and checkout relevant version
RUN git clone https://www.github.com/coin-or/Ipopt.git $IPOPT_WITH_VERS
WORKDIR $IPOPT_WITH_VERS
RUN git checkout releases/$IPOPT_VERS

# Make a build directory and configure
RUN mkdir -p build
WORKDIR build

RUN ../configure --with-hsl --with-lapack --with-lapack-lflags="-lf77blas -lcblas -latlas -llapack" && \
  make && \
  make test && \
  make install

RUN sudo apt install -yq libomp-dev

#=============================================================
# CUDA BEGIN
#=============================================================

ENV NVARCH x86_64

ENV NVIDIA_REQUIRE_CUDA "cuda>=12.1 brand=tesla,driver>=450,driver<451 brand=tesla,driver>=470,driver<471 brand=unknown,driver>=470,driver<471 brand=nvidia,driver>=470,driver<471 brand=nvidiartx,driver>=470,driver<471 brand=geforce,driver>=470,driver<471 brand=geforcertx,driver>=470,driver<471 brand=quadro,driver>=470,driver<471 brand=quadrortx,driver>=470,driver<471 brand=titan,driver>=470,driver<471 brand=titanrtx,driver>=470,driver<471 brand=tesla,driver>=510,driver<511 brand=unknown,driver>=510,driver<511 brand=nvidia,driver>=510,driver<511 brand=nvidiartx,driver>=510,driver<511 brand=geforce,driver>=510,driver<511 brand=geforcertx,driver>=510,driver<511 brand=quadro,driver>=510,driver<511 brand=quadrortx,driver>=510,driver<511 brand=titan,driver>=510,driver<511 brand=titanrtx,driver>=510,driver<511 brand=tesla,driver>=515,driver<516 brand=unknown,driver>=515,driver<516 brand=nvidia,driver>=515,driver<516 brand=nvidiartx,driver>=515,driver<516 brand=geforce,driver>=515,driver<516 brand=geforcertx,driver>=515,driver<516 brand=quadro,driver>=515,driver<516 brand=quadrortx,driver>=515,driver<516 brand=titan,driver>=515,driver<516 brand=titanrtx,driver>=515,driver<516 brand=tesla,driver>=525,driver<526 brand=unknown,driver>=525,driver<526 brand=nvidia,driver>=525,driver<526 brand=nvidiartx,driver>=525,driver<526 brand=geforce,driver>=525,driver<526 brand=geforcertx,driver>=525,driver<526 brand=quadro,driver>=525,driver<526 brand=quadrortx,driver>=525,driver<526 brand=titan,driver>=525,driver<526 brand=titanrtx,driver>=525,driver<526"
ENV NV_CUDA_CUDART_VERSION 12.1.55-1
ENV NV_CUDA_COMPAT_PACKAGE cuda-compat-12-1

ENV NVARCH x86_64
ENV NVIDIA_REQUIRE_CUDA "cuda>=12.1"
ENV NV_CUDA_CUDART_VERSION 12.1.55-1

RUN apt-get update && apt-get install -y --no-install-recommends \
    gnupg2 curl ca-certificates && \
    curl -fsSL https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/${NVARCH}/3bf863cc.pub | apt-key add - && \
    echo "deb https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/${NVARCH} /" > /etc/apt/sources.list.d/cuda.list #&& \
    #apt-get purge --autoremove -y curl # \
    # && rm -rf /var/lib/apt/lists/*

ENV CUDA_VERSION 12.1.0

## For libraries in the cuda-compat-* package: https://docs.nvidia.com/cuda/eula/index.html#attachment-a
#RUN apt-get update && apt-get install -y --no-install-recommends \
#    cuda-cudart-12-1=${NV_CUDA_CUDART_VERSION} \
#    ${NV_CUDA_COMPAT_PACKAGE} # \
#    #&& rm -rf /var/lib/apt/lists/*

# Required for nvidia-docker v1
RUN echo "/usr/local/nvidia/lib" >> /etc/ld.so.conf.d/nvidia.conf \
    && echo "/usr/local/nvidia/lib64" >> /etc/ld.so.conf.d/nvidia.conf

ENV PATH /usr/local/nvidia/bin:/usr/local/cuda/bin:${PATH}
ENV LD_LIBRARY_PATH /usr/local/nvidia/lib:/usr/local/nvidia/lib64

# nvidia-container-runtime
ENV NVIDIA_VISIBLE_DEVICES all
ENV NVIDIA_DRIVER_CAPABILITIES compute,utility

ENV NV_CUDA_LIB_VERSION "12.1.0-1"

ENV NV_CUDA_CUDART_DEV_VERSION 12.1.55-1
ENV NV_NVML_DEV_VERSION 12.1.55-1
ENV NV_LIBCUSPARSE_DEV_VERSION 12.0.2.55-1
ENV NV_LIBNPP_DEV_VERSION 12.0.2.50-1
ENV NV_LIBNPP_DEV_PACKAGE libnpp-dev-12-1=${NV_LIBNPP_DEV_VERSION}

ENV NV_LIBCUBLAS_DEV_VERSION 12.1.0.26-1
ENV NV_LIBCUBLAS_DEV_PACKAGE_NAME libcublas-dev-12-1
ENV NV_LIBCUBLAS_DEV_PACKAGE ${NV_LIBCUBLAS_DEV_PACKAGE_NAME}=${NV_LIBCUBLAS_DEV_VERSION}

ENV NV_NVPROF_VERSION 12.1.55-1
ENV NV_NVPROF_DEV_PACKAGE cuda-nvprof-12-1=${NV_NVPROF_VERSION}

ENV NV_LIBNCCL_DEV_PACKAGE_NAME libnccl-dev
ENV NV_LIBNCCL_DEV_PACKAGE_VERSION 2.17.1-1
ENV NCCL_VERSION 2.17.1-1
ENV NV_LIBNCCL_DEV_PACKAGE ${NV_LIBNCCL_DEV_PACKAGE_NAME}=${NV_LIBNCCL_DEV_PACKAGE_VERSION}+cuda12.1

ENV NV_CUDA_CUDART_DEV_VERSION 12.1.55-1
ENV NV_NVML_DEV_VERSION 12.1.55-1
ENV NV_LIBCUSPARSE_DEV_VERSION 12.0.2.55-1
ENV NV_LIBNPP_DEV_VERSION 12.0.2.50-1
ENV NV_LIBNPP_DEV_PACKAGE libnpp-dev-12-1=${NV_LIBNPP_DEV_VERSION}

ENV NV_LIBCUBLAS_DEV_PACKAGE_NAME libcublas-dev-12-1
ENV NV_LIBCUBLAS_DEV_VERSION 12.1.0.26-1
ENV NV_LIBCUBLAS_DEV_PACKAGE ${NV_LIBCUBLAS_DEV_PACKAGE_NAME}=${NV_LIBCUBLAS_DEV_VERSION}

ENV NV_LIBNCCL_DEV_PACKAGE_NAME libnccl-dev
ENV NV_LIBNCCL_DEV_PACKAGE_VERSION 2.17.1-1
ENV NCCL_VERSION 2.17.1-1
ENV NV_LIBNCCL_DEV_PACKAGE ${NV_LIBNCCL_DEV_PACKAGE_NAME}=${NV_LIBNCCL_DEV_PACKAGE_VERSION}+cuda12.1


ARG TARGETARCH

RUN apt-get update 
RUN apt-get install libnccl2=2.17.1-1+cuda12.1
RUN apt-get install -y --no-install-recommends \
    libtinfo5 libncursesw5 \
    cuda-cudart-dev-12-1=${NV_CUDA_CUDART_DEV_VERSION} \
    cuda-command-line-tools-12-1=${NV_CUDA_LIB_VERSION} \
    cuda-minimal-build-12-1=${NV_CUDA_LIB_VERSION} \
    cuda-libraries-dev-12-1=${NV_CUDA_LIB_VERSION} \
    cuda-nvml-dev-12-1=${NV_NVML_DEV_VERSION} \
    ${NV_NVPROF_DEV_PACKAGE} \
    ${NV_LIBNPP_DEV_PACKAGE} \
    libcusparse-dev-12-1=${NV_LIBCUSPARSE_DEV_VERSION} \
    ${NV_LIBCUBLAS_DEV_PACKAGE} \
    ${NV_LIBNCCL_DEV_PACKAGE} 

#\
#&& rm -rf /var/lib/apt/lists/*

# Keep apt from auto upgrading the cublas and nccl packages. See https://gitlab.com/nvidia/container-images/cuda/-/issues/88

RUN apt-mark hold ${NV_LIBCUBLAS_DEV_PACKAGE_NAME} ${NV_LIBNCCL_DEV_PACKAGE_NAME}
ENV LIBRARY_PATH /usr/local/cuda/lib64/stubs



#=============================================================
# CUDA END
#=============================================================










RUN apt-get install -yq ssh rsync
RUN apt-get install -yq libc++-12-dev libc++abi-12-dev libclang-12-dev libomp5-12 libomp-12-dev llvm-12-dev llvm-12-tools
RUN wget https://apt.llvm.org/llvm.sh && chmod +x llvm.sh && ./llvm.sh 15 all

RUN apt-get -yq update
#RUN apt-get install -yq gnuplot-x11
#RUN apt-get install -yq gnuplot

RUN apt install -yq qtbase5-dev qtbase5-dev-tools libqt5svg5-dev qttools5-dev-tools
RUN apt install -yq strace

WORKDIR $TOP_LEVEL_SRC_DIR
RUN wget https://jaist.dl.sourceforge.net/project/gnuplot/gnuplot/5.4.1/gnuplot-5.4.1.tar.gz
RUN tar -xvzf gnuplot-5.4.1.tar.gz
WORKDIR $TOP_LEVEL_SRC_DIR/gnuplot-5.4.1
RUN ./configure --with-qt=qt5 && \
    make -j && \
    make install

#RUN apt-get install -yq libc++-14-dev libc++abi-14-dev libclang-14-dev libomp5-14 libomp-14-dev llvm-14-dev llvm-14-tools
#RUN apt-get install -yq clang++-14
#RUN wget https://apt.llvm.org/llvm.sh && chmod +x llvm.sh && ./llvm.sh 15 all
#RUn apt search format | grep clang
WORKDIR $TOP_LEVEL_SRC_DIR/iwyu-dir
RUN git clone https://github.com/include-what-you-use/include-what-you-use.git && \
    cd include-what-you-use && \
    git checkout clang_15
WORKDIR $TOP_LEVEL_SRC_DIR/iwyu-dir/build
RUN cmake -G "Unix Makefiles" -DCMAKE_PREFIX_PATH=/usr/lib/llvm-15 ../include-what-you-use
RUN make -j install

WORKDIR $TOP_LEVEL_SRC_DIR/mold-linker
RUN git clone https://github.com/rui314/mold.git && \
  mkdir mold/build && \
  cd mold/build && \
  git checkout v1.7.1 && \
  ../install-build-deps.sh && \
  cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=clang++-15 .. && \
  cmake --build . -j $(nproc) && \
  cmake --install .

WORKDIR $TOP_LEVEL_SRC_DIR

#RUN sh -c "$(wget -O- https://github.com/deluan/zsh-in-docker/releases/download/v1.1.2/zsh-in-docker.sh)"

ENV PATH /usr/local/nvidia/bin:/usr/local/cuda-11.7/bin:${PATH}
#ENV LD_LIBRARY_PATH /usr/local/nvidia/lib:/usr/local/nvidia/lib64

#=============================================================
# CUDA END
#=============================================================

RUN wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip && \
    unzip eigen-3.4.0.zip
RUN rm -rf /usr/include/eigen3 && \
    cd eigen-3.4.0 && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=/usr .. && \
    make install

# ----------------
# Install valgrind
# ----------------

RUN apt-get install -yq --allow-downgrades \
  libc6=2.31-0ubuntu9.9 \
  libc6-dbg=2.31-0ubuntu9.9 \
  libc-dev-bin=2.31-0ubuntu9.9 \
  libc6-dev=2.31-0ubuntu9.9 \
  libc6-dev-i386=2.31-0ubuntu9.9 \
  libc6-dev-x32=2.31-0ubuntu9.9 \
  libc6-i386=2.31-0ubuntu9.9 \
  libc6-pic=2.31-0ubuntu9.9 \
  libc6-prof=2.31-0ubuntu9.9 \
  libc6-x32=2.31-0ubuntu9.9

RUN apt-get install -yq linux-tools-5.15.0-56-generic linux-cloud-tools-5.15.0-56-generic linux-tools-generic linux-tools-common
RUN apt-get install -yq iputils-ping
#RUN apt-get install -yq linux-tools-generic linux-cloud-tools-generic


# Download the valgrind package manually, install with no-deps, and remove installation file
RUN wget http://archive.ubuntu.com/ubuntu/pool/main/v/valgrind/valgrind_3.15.0-1ubuntu9.1_amd64.deb && \
 dpkg -i valgrind_3.15.0-1ubuntu9.1_amd64.deb && \
 rm -f valgrind_3.15.0-1ubuntu9.1_amd64.deb
#RUN wget http://archive.ubuntu.com/ubuntu/pool/main/v/valgrind/valgrind_3.15.0-1ubuntu9.1_amd64.deb && \
# dpkg --ignore-depends=libc6-dbg -i valgrind_3.15.0-1ubuntu9.1_amd64.deb && \
# rm -f valgrind_3.15.0-1ubuntu9.1_amd64.deb

ENV PATH="${PATH}:/usr/local/cuda-11.7/bin"

WORKDIR /installs
RUN git clone https://github.com/ninja-build/ninja.git && \
    cd ninja && \
    git checkout release && \
    cmake -Bbuild-cmake && \
    cmake --build build-cmake
RUN cd ninja && \
    cd build-cmake && \
    make install

WORKDIR /installs/bazel-installs
RUN wget https://github.com/bazelbuild/bazel/releases/download/6.2.0/bazel-6.2.0-linux-x86_64
RUN chmod a+x bazel-6.2.0-linux-x86_64

WORKDIR /simulator
#RUN /bin/bash -c "cd /installs/cpp-opt/src && source /opt/ros/noetic/setup.bash && catkin build"

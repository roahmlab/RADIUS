CONTAINER_NAME="RADIUS_simulator"

CONTAINER_ID=`docker ps -aqf "name=^/${CONTAINER_NAME}$"`
if [ -z "${CONTAINER_ID}" ]; then
  echo "Creating new container ${CONTAINER_NAME}."
  cp /etc/passwd "$(pwd)/.etc_passwd"
  getent passwd $(whoami) >> "$(pwd)/.etc_passwd"

  cp /etc/group "$(pwd)/.etc_group"
  getent group "$(whoami)" >> "$(pwd)/.etc_group"
    #-u 0:0                   \

  xhost +
  docker run -it --rm                     \
    --privileged \
    --cap-add SYS_ADMIN \
    -e DISPLAY=$DISPLAY                    \
    -v /tmp/.X11-unix:/tmp/.X11-unix:ro    \
    --cap-add sys_ptrace                   \
    -p127.0.0.1:2222:22                    \
    --name=${CONTAINER_NAME}               \
    --gpus=all                             \
    --shm-size=512M                        \
    --network="host"                       \
    -u $(id -u):$(id -g)                   \
    -ti                                    \
    -v "$(pwd)/.etc_passwd":/etc/passwd:ro   \
    -v "$(pwd)/.etc_group":/etc/group:ro     \
    -v "$(pwd)/docker-home":"/home/$(id -nu)"  \
    -v /dev/shm:/dev/shm                  \
    -v /etc/timezone:/etc/timezone:ro     \
    -v /etc/localtime:/etc/localtime:ro   \
    -v $(pwd)/../util:/data:rw            \
    -v $(pwd)/simulator:/simulator        \
    -v "$(pwd)/..":"/home/$(id -nu)/radius-repo":rw \
    -w /simulator                         \
  roahm/radius_sim /usr/bin/bash
else
  if [ -z `docker ps -qf "name=^/${CONTAINER_NAME}$"` ]; then
    echo "Found ${CONTAINER_NAME} docker container: ${CONTAINER_ID}."
    xhost +local:${CONTAINER_ID}
    echo "Starting and attaching to ${CONTAINER_NAME} container..."
    docker start ${CONTAINER_ID}
    docker attach ${CONTAINER_ID}
  else
    echo "Found running ${CONTAINER_NAME} container, attaching bash..."
    docker exec -it ${CONTAINER_ID} bash
  fi
fi

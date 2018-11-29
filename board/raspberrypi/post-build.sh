#!/bin/sh

set -u
set -e

# Add a console on tty1
if [ -e ${TARGET_DIR}/etc/inittab ]; then
    grep -qE '^tty1::' ${TARGET_DIR}/etc/inittab || \
	sed -i '/GENERIC_SERIAL/a\
tty1::respawn:/sbin/getty -L  tty1 0 vt100 # HDMI console' ${TARGET_DIR}/etc/inittab
fi

  
cp $BASE_DIR/../custom-scripts/S41network-config $BASE_DIR/target/etc/init.d
chmod +x $BASE_DIR/target/etc/init.d/S41network-config

cp $BASE_DIR/../custom-scripts/pyServer.py $BASE_DIR/target/etc/init.d
chmod +x $BASE_DIR/target/etc/init.d/pyServer.py

cp $BASE_DIR/../custom-scripts/deadline $BASE_DIR/target/etc/init.d
chmod +x $BASE_DIR/target/etc/init.d/deadline

cp $BASE_DIR/../custom-scripts/FalcOpt_dU/mainDeadline $BASE_DIR/target/etc/init.d
chmod +x $BASE_DIR/target/etc/init.d/mainDeadline

cp $BASE_DIR/../custom-scripts/FalcOpt_dU/sim -r $BASE_DIR/target/etc/init.d
chmod +x $BASE_DIR/target/etc/init.d/sim

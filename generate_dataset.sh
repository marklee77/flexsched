#!/bin/bash

if [ $# != 1 ]; then
  echo "usage: $0 <basedir name>"
  exit 1
fi

BASEDIR=$1

NUMINSTANCES=100
NUMRIGID=2
NUMFLUID=2
NUMSERVERS=64
NUMSERVICES=100
RIGIDSLACK="0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9"
RIGIDSLACK="0.5"
RIGIDCV=".5"
FLUIDMEAN=".5"
FLUIDCV=".5"
SLAPROB="0.0 0.5 1.0"
SLAVALUE=".5"

\rm -rf $BASEDIR
echo "./genprobset.sh $BASEDIR \"$NUMINSTANCES\" \"$NUMRIGID\" \"$NUMFLUID\" \"$NUMSERVERS\" \"$NUMSERVICES\" \"$RIGIDSLACK\" \"$RIGIDCV\" \"$FLUIDMEAN\" \"$FLUIDCV\" \"$SLAPROB\" \"$SLAVALUE\""
./genprobset.sh $BASEDIR "$NUMINSTANCES" "$NUMRIGID" "$NUMFLUID" "$NUMSERVERS" "$NUMSERVICES" "$RIGIDSLACK" "$RIGIDCV" "$FLUIDMEAN" "$FLUIDCV" "$SLAPROB" "$SLAVALUE"


echo "Number of files (each file = one problem instance): $NUMINSTANCES">> $BASEDIR/README
echo "">> $BASEDIR/README
echo "File Names:">> $BASEDIR/README
echo "Each file is named problem-a-b-c-d-e-f-g-h-i-j-k-l.txt with:">> $BASEDIR/README
echo "a = number of rigid resources ($NUMRIGID)">> $BASEDIR/README
echo "b = number of fluid resources ($NUMFLUID)">> $BASEDIR/README
echo "c = number of servers ($NUMSERVERS)">> $BASEDIR/README
echo "d = number of services ($NUMSERVICES)">> $BASEDIR/README
echo "e = rigid slack ($RIGIDSLACK)">> $BASEDIR/README
echo "f = rigid c.v. ($RIGIDCV)">> $BASEDIR/README
echo "g = mean fluid need ($FLUIDMEAN)">> $BASEDIR/README
echo "h = fluid c.v. ($FLUIDCV)">> $BASEDIR/README
echo "i = sla prob  ($SLAPROB)">> $BASEDIR/README
echo "j = sla value  ($SLAVALUE)">> $BASEDIR/README
echo "k = rng seed used">> $BASEDIR/README
echo "l = instance number">> $BASEDIR/README
echo "">> $BASEDIR/README
echo "">> $BASEDIR/README
echo "">> $BASEDIR/README
echo "File Content:">> $BASEDIR/README
echo "<number of rigid needs>">> $BASEDIR/README
echo "<number of fluid needs>">> $BASEDIR/README
echo "<number of servers>">> $BASEDIR/README
echo "<number of services>">> $BASEDIR/README
echo "<sla> <rigid 1 of svc 1> <rigid 2 of svc 1> ... <fluid 1 of svc 1> <fluid 2 of svc 1>...">> $BASEDIR/README
echo "<sla> <rigid 1 of svc 2> <rigid 2 of svc 2> ... <fluid 1 of svc 2> <fluid 2 of svc 2>...">> $BASEDIR/README
echo "....">> $BASEDIR/README

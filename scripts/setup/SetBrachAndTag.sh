export CATVERSION=v7-6-3
### If there is a small bug/new code then new subtag is made
tag_numerator='.4'
if [[ '-d' == branch ]];
    then
    export CATTAG=
else
    export CATTAG=v7-6-3.4
fi


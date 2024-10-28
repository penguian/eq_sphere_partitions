#!/bin/bash
pathf=$1
f=$(basename ${pathf/.m/})
F=${f^^}
patht="eq_test/test_${f}.m"
t=$(basename ${patht/.m/})
T=${t^^}
Revision=${2:-1.12}
Date=$(date -I)

echo $t
cat << END > $patht
function $t
%$T Test the $f function

% Copyright 2024 Paul Leopardi.
% \$Revision \$ $Revision \$ \$Date $Date \$
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help $f
END
grep "^% >> " $pathf >> $patht
sed -i 's/% >> //' $patht

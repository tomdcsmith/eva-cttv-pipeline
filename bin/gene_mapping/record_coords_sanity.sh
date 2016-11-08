#!/usr/bin/env bash

## $1 is variant summary file
## $2 is coordinate file (non compressed)

# compare RCVs
diff <(zcat $1 | grep -i pathogenic | grep -i grch38 | cut -f9 | tr ';' '\n' | sort | uniq) <(cut -f7 $2 | sort | uniq)

# compare rs ids
diff <(zcat $1 | grep -i pathogenic | grep -i grch38 | cut -f7 | tr ';' '\n' | sort | uniq) <(cut -f6 $2 | sort | uniq | sed 's/rs//')

# compare nsv ids
diff <(zcat $1 | grep -i pathogenic | grep -i grch38 | cut -f8 | tr ';' '\n' | sort | uniq) <(cut -f9 $2 | sort | uniq)




#!/usr/bin/env python3
"""Read-only diagnosis of accepted-event ordering; does not approve closure."""
import collections
import json
from pathlib import Path
import sys


def keyed(path):
    result=[]
    with Path(path).open() as stream:
        for line in stream:
            record=json.loads(line)
            key=tuple(tuple(leg[:3])+tuple(round(x,8) for x in leg[3:]) for leg in record["hard"])
            result.append((key,record["sequence"]))
    return result


def main():
    a,b=(keyed(path) for path in sys.argv[1:3])
    ca,cb=(collections.Counter(key for key,_ in sample) for sample in (a,b))
    lookup={key:seq for key,seq in a}
    paired=[(lookup[key],seq) for key,seq in b if key in lookup]
    first_mismatch=next(((i+1,ka,kb) for i,((ka,_),(kb,_)) in enumerate(zip(a,b)) if ka!=kb),None)
    print(json.dumps({"reference_count":len(a),"candidate_count":len(b),
        "reference_unique":len(ca),"candidate_unique":len(cb),
        "reference_missing":sum((ca-cb).values()),"candidate_extra":sum((cb-ca).values()),
        "candidate_keys_not_in_reference":len(set(cb)-set(ca)),
        "first_order_mismatch":first_mismatch,
        "offset_min":min((x-y for x,y in paired),default=None),
        "offset_max":max((x-y for x,y in paired),default=None),
        "most_common_offsets":collections.Counter(x-y for x,y in paired).most_common(8)},indent=2))


if __name__ == "__main__":
    main()

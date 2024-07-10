"""Remove failed WorkChainNode from a group"""
import sys

from aiida import orm
from aiida import load_profile

def main():
    load_profile()

    # only support group name
    group = orm.load_group(label=sys.argv[1])

    nodes_to_delete = []
    for node in group.nodes:
        if not isinstance(node, orm.WorkChainNode):
            continue

        if not node.is_finished_ok:
            nodes_to_delete.append(node)

    if len(sys.argv) > 2 and sys.argv[2] == '--dry':
        print(f"Going to remove {[i.pk for i in nodes_to_delete]}.") 
    else:
        group.remove_nodes(nodes_to_delete)
        print(f"{[i.pk for i in nodes_to_delete]} removed from group.") 

if __name__ == '__main__':
    main()

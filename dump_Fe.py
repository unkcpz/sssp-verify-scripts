import json
import os

from aiida import orm
from aiida.orm import load_node

def flaten_output(attr_dict):
    """flaten output dict node"""
    for key, value in attr_dict.items():
        if isinstance(value, orm.Dict):
            attr_dict[key] = value.get_dict()
        else:
            flaten_output(value)

def dump_output(element, pks, folder):
    """dump the verification result
    if fname is None, create the name from label
    """
    res = {}
    for pk in pks:
        node = load_node(pk)
        label = node.extras.get('label')
        
        output_dict = {}
        outputs = node.outputs._construct_attribute_dict(False)
        flaten_output(outputs)
        
        output_dict.update(outputs)
        
        res[label] = output_dict
    
    abs_path = os.path.join(folder, f'{element}.json')

    with open(abs_path, 'w') as fh:
        json.dump(dict(res), fh, indent=2)
    
if __name__ == '__main__':
    # TODO: add a function to search the node and dry_run it with display the node info
    # wait for determination. 
    pk_list = [112579, 112588, 112597, 112606, 112615, 112624, 112633, 112642, 112651, 112660]
    dump_output(element='Fe', pks=pk_list, folder='./result_json/')
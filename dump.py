import json
import os

import click
from aiida import orm
from aiida.orm import load_node

def flaten_output(attr_dict):
    """flaten output dict node"""
    for key, value in attr_dict.items():
        if isinstance(value, orm.Dict):
            attr_dict[key] = value.get_dict()
        else:
            flaten_output(value)

def dump_output(element, pks, folder, append):
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

    # with open(abs_path, 'r') as fh:
    #     data = json.load(fh)
    #     for k, v in res.items():
    #         if k in data and override:
    #             data[k] = v
                
    #         if k not in data:
    #             data[k] = v
                
    with open(abs_path, 'w') as fh:
        json.dump(dict(res), fh, indent=2)
        
@click.command()
@click.option('--element', help='element to dump.')
@click.option('-f', '--folder', default='./result_json/', help='the folder to store the json file.')
@click.option('-p', '--append/--no-append', default=False, help='if the pseudo info exist append it.')
@click.argument('nodes', nargs=-1)
def run(element, folder, append, nodes):
    # pk_list = [62683, 62693, 62709, 62727, 62814, 62895, 62968, 63045]
    dump_output(element=element, pks=nodes, folder=folder, append=append)
    
if __name__ == '__main__':
    # TODO: add a function to search the node and dry_run it with display the node info
    # wait for determination. 
    run()
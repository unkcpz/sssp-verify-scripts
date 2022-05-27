"""Giving the node invalid cache if it is a calcjob or 
of its descendant calcjob if it is a workchain."""

import click

import aiida
from aiida import orm

aiida.load_profile()

@click.command()
@click.option('profile', '-p', help='profile')
@click.option('--find-all', is_flag=True, default=True, 
              help='also find same nodes (cached or source) and invalid_cache.')
@click.argument('nodes', type=int, nargs=-1)
def run(profile, find_all, nodes):
    _profile = aiida.load_profile(profile)
    click.echo(f'Profile: {_profile.name}')
    
    nodes_to_process = []
    for n in nodes:
        n = orm.load_node(n)
        if isinstance(n, orm.CalcJobNode):
            nodes_to_process.append(n)
        elif isinstance(n, orm.WorkChainNode):
            for called_descendant in n.called_descendants:
                if isinstance(called_descendant, orm.CalcJobNode):
                    nodes_to_process.append(called_descendant)
        else:
            raise TypeError(f"{type(n)} not supported.")
        
    for n in nodes_to_process:
        if find_all:
            same_nodes = n.get_all_same_nodes
        else:
            same_nodes = [n]
        for nn in same_nodes:
            try:
                click.echo(f'cleanning hash of node pk={nn.pk}')
                nn.delete_extra('_aiida_hash')
            except:
                click.echo(f'{nn} do not have `_aiida_hash` field, might be already cleaned.')

if __name__ == '__main__':
    run()
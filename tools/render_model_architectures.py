import argparse
from pathlib import Path
from diagrams import Diagram, Node, Edge
p=argparse.ArgumentParser(); p.add_argument('--outdir',required=True); a=p.parse_args()
out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
for model in ('m2a','m2b','m3a','m3b'):
    with Diagram(model.upper()+' | Candidate-ranking model',filename=str(out/model),outformat=['svg','png'],show=False,direction='TB',graph_attr={'bgcolor':'white','fontname':'Helvetica','fontsize':'21','labelloc':'t','pad':'0.45','ranksep':'0.7','nodesep':'0.5'},edge_attr={'color':'#7b8799','arrowsize':'0.8'}):
        def box(label,color='#168fa5'):
            return Node(label,shape='box',style='rounded,filled',fillcolor='white',color=color,penwidth='1.8',labelloc='c',fontname='Helvetica',fontsize='12',fixedsize='false',margin='0.22,0.17',width='2.8',height='0.9')
        x=box('Prepared candidate features\nX and case boundaries')
        priors=box('Coefficient priors\nNormal population coefficients','#8250c4')
        if model.startswith('m3'):
            groups=box('Partial pooling\nbeta + tau_beta * z_beta','#8250c4')
            priors >> groups
            linear=box('Group-specific linear score\nX times group coefficients')
            groups >> linear
        else:
            linear=box('Global linear score\nX times beta')
            priors >> linear
        x >> linear
        if model.endswith('b'):
            u=box('HPO and network transforms\nPseudo-observations + availability mask')
            if model=='m3b':
                cp=box('Hierarchical dependence\nGroup correlation and copula coefficient','#8250c4')
            else:
                cp=box('Global dependence\nrho = tanh(rho_raw); beta_cop','#8250c4')
            cop=box('Gaussian-copula contribution\nMasked when unavailable')
            u >> cop; cp >> cop
            eta=box('Combined candidate score\nLinear + copula contribution')
            linear >> eta; cop >> eta
        else: eta=linear
        likelihood=box('Training likelihood\nCategorical-logit within each case','#159570')
        prediction=box('New-case prediction\nWithin-case softmax gives p_new','#159570')
        eta >> likelihood
        eta >> prediction

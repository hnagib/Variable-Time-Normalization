import statsmodels.formula.api as smf
import pandas as pd
import sys
from bokeh.plotting import figure, show
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.io import output_notebook
from bokeh.palettes import Dark2_5 as palette
from bokeh.layouts import row
import itertools 
import numpy as np

class VariableTimeNormalization:
    '''
    Variable time normalization analysis
    '''
    def __init__(self, df0, exp_pairs, rxn_orders, product_name, kobs_line_format='y ~ x-1'):
        
        self.df0 = df0
        self.exp_pairs = exp_pairs
        self.rxn_orders = rxn_orders
        self.product_name = product_name
        self.kobs_line_format = kobs_line_format

    @property
    def var_conc_comp(self):
        var_conc_comp=[]
        for key in self.exp_pairs.keys():
            if self.exp_pairs[key]['type'] == 'variable':
                var_conc_comp.append(key)
        return var_conc_comp
    
    @property
    def const_conc_comp(self):
        const_conc_comp=[]
        for key in self.exp_pairs.keys():
            if self.exp_pairs[key]['type'] == 'constant':
                const_conc_comp.append(key)
        return const_conc_comp
    
    @property
    def kobs_normalizer(self):
        return '∑'+str.join('',[f'[{col}]^{self.rxn_orders[col]}' for col in self.var_conc_comp + self.const_conc_comp])+'∆t'
    
    @property
    def df_vtn(self):
        df = self.df0
        for key in df.keys():

            kobs_cols = []
            for col in df[key].columns:

                # Add time lagged columns
                if col in self.var_conc_comp + ['t']:
                    df[key][f'{col}-1'] = df[key][col].shift(1)

                # Add three columns for each reactant: [A]^α, ([A]^α)*∆t and ∑([A]^α)*∆t        
                if col in self.var_conc_comp:
                    kobs_cols.append(f'[{col}]^{self.rxn_orders[col]}')
                    df[key][f'[{col}]^{self.rxn_orders[col]}'] = 0.5*((df[key][f'{col}'] + df[key][f'{col}-1']).pow(self.rxn_orders[col]))
                    df[key][f'[{col}]^{self.rxn_orders[col]}∆t'] = 0.5*((df[key][f'{col}'] + df[key][f'{col}-1']).pow(self.rxn_orders[col]))*(df[key]['t'] - df[key]['t-1'])
                    df[key][f'∑[{col}]^{self.rxn_orders[col]}∆t'] = df[key][f'[{col}]^{self.rxn_orders[col]}∆t'].cumsum()

                # Add two columns for each catalyst: [cat]^γ, t*[cat]^γ
                for cat_name in self.const_conc_comp:
                    kobs_cols.append(f'[{cat_name}]^{self.rxn_orders[cat_name]}')
                    df[key][f'[{cat_name}]^{self.rxn_orders[cat_name]}'] = df[key][cat_name].pow(self.rxn_orders[cat_name])
                    df[key][f't[{cat_name}]^{self.rxn_orders[cat_name]}'] = df[key]['t']*df[key][cat_name].pow(self.rxn_orders[cat_name])

            # Add time delta column: ∆t
            kobs_cols.append('∆t')
            df[key]['∆t'] = df[key]['t'] - df[key]['t-1']

            # Add column used for computing kobs: ∑[A]^α * [B]^β * [cat]^γ * ∆t
            df[key][self.kobs_normalizer] = df[key][list(set(kobs_cols))].product(axis=1, skipna=False).cumsum()
        return df
    
    @property
    def kobs_cols(self):
        df = self.df0
        for key in df.keys():

            kobs_cols = []
            for col in df[key].columns:
                # Add three columns for each reactant: [A]^α, ([A]^α)*∆t and ∑([A]^α)*∆t        
                if col in self.var_conc_comp:
                    kobs_cols.append(f'[{col}]^{self.rxn_orders[col]}')
                # Add two columns for each catalyst: [cat]^γ, t*[cat]^γ
                for cat_name in self.const_conc_comp:
                    kobs_cols.append(f'[{cat_name}]^{self.rxn_orders[cat_name]}')
            # Add time delta column: ∆t
            kobs_cols.append('∆t')
        return list(set(kobs_cols))
    
    def plot_vtn(self, comp):
        plot = figure(width=800,height=350)
        
        if comp in self.var_conc_comp:
            x_var="∑[{comp}]^{rxn_order}∆t".format(comp=comp, rxn_order=self.rxn_orders[comp])
        
        if comp in self.const_conc_comp:
            x_var="t[{comp}]^{rxn_order}".format(comp=comp, rxn_order=self.rxn_orders[comp])
        
        for exp, color in zip(self.exp_pairs[comp]['exps'], self.exp_pairs[comp]['colors']):
            plot.circle(
                x=x_var, 
                y=self.product_name, 
                color=color, 
                source=ColumnDataSource(self.df_vtn[exp]),
                legend=f'Experiment-{exp}',
                size=7.5,
                alpha=0.5
            )
        
        x = pd.concat([
            self.df_vtn[exp][[self.product_name,x_var]] 
            for exp in self.exp_pairs[comp]['exps']]
        )
        plot.line(x=x_var, y=self.product_name, source=x.sort_values(x_var), legend='line', color='grey')
        
        plot.add_tools(HoverTool(tooltips=[(i, f'@{i}'+'{0.000}') for i in [self.product_name]]))
        plot.xaxis.axis_label = x_var
        plot.yaxis.axis_label = f"[{self.product_name}]"
        plot.legend.location = "top_left"
        plot.legend.click_policy="hide"
        show(plot)
    
    @property
    def df_kobs_reg(self):
        df = pd.concat([self.df_vtn[f'{exp}'][[self.product_name, self.kobs_normalizer]] for exp in self.df0.keys()]).dropna()
        df['y'] = df[self.product_name]
        df['x'] = df[self.kobs_normalizer]
        return df
    
    @property
    def kobs_line(self):
        return smf.ols(self.kobs_line_format, data=self.df_kobs_reg[['x','y']]).fit()
    
    @property
    def kobs(self):
        return self.kobs_line.params['x']
    
    def plot_kobs(self):
        plot = figure(width=800,height=350, title=f'kobs = {round(self.kobs,3)}')
        x_var = self.kobs_normalizer
        for key, color in zip(self.df0.keys(), itertools.cycle(palette)):
            plot.circle(
                x=x_var, 
                y=self.product_name, 
                color=color, 
                source=ColumnDataSource(self.df_vtn[key]),
                legend=f'Experiment-{key}',
                size=7.5,
                alpha=0.5
            )
        plot.line(
            x='X',
            y='Yp',
            source=pd.DataFrame({'X':self.df_kobs_reg['x'], 'Yp':self.kobs_line.predict(self.df_kobs_reg['x'])}),
            color='grey',
            legend=f'ols: {self.kobs_line_format}'
        )
        plot.xaxis.axis_label = x_var
        plot.yaxis.axis_label = f"[{self.product_name}]"
        
        plot.legend.location = "top_left"
        plot.legend.click_policy="hide"
        show(plot)
        
def tv(x0, exp, comp):
    x0 = x0[0]
    exp.rxn_orders[comp]=x0
    
    if comp in exp.var_conc_comp:
        x_var="∑[{comp}]^{rxn_order}∆t"

    if comp in exp.const_conc_comp:
        x_var="t[{comp}]^{rxn_order}"
        
    x = pd.concat([
        exp.df_vtn[i][[exp.product_name,x_var.format(comp=comp, rxn_order=x0)]] 
        for i in exp.exp_pairs[comp]['exps']]
    )
    
    x.columns=[exp.product_name,'t']
    x=x.sort_values('t').dropna()
    x[f'{exp.product_name}_diff'] = x[exp.product_name].diff(1)
    x['t'] = x['t']/x['t'].max()
    x['t_diff'] = x['t'].diff(1)
    x['tv'] = np.sqrt(x[f'{exp.product_name}_diff'].pow(2)+x['t_diff'].pow(2))
    return x['tv'].sum()
        
def plot_tv(exp, comps, x0_arange):
    plots=[]
    for comp in comps:
        p = []
        for i in np.arange(x0_arange[0],x0_arange[1], x0_arange[2]):
            p.append([i, tv([i], exp, comp)])
        df = pd.DataFrame(p)
        plot = figure(width=325,height=300, title=f'{comp}')
        x=f'Order of [{comp}]'
        y='TV()'
        df.columns=[x,y]
        plot.line(x, y, source=ColumnDataSource(df), color='grey')    
        plot.xaxis.axis_label = x
        plot.yaxis.axis_label = y
        #plot.add_tools(HoverTool(tooltips=[(i, f'@{i}'+'{0.000}') for i in [x, y]]))
        plots.append(plot)
    
    show(row(plots))
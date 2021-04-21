library(DiagrammeR)
library(DiagrammeRsvg)
t<-DiagrammeR::grViz("digraph {

graph [layout = dot, rankdir = TB]

# define the global styles of the nodes. We can override these in box if we wish
node [shape = rectangle, style = filled, fillcolor = Linen]

data1 [label = 'Baseline Data', shape = folder, fillcolor = Turquoise]
a1 [label =  'Households \n n = 2994\nIndividuals \n n = 5313', fillcolor = Turquoise]

a6 [label = 'Individuals \n n = 869']
a7 [label = 'Individuals \n n = 4444']
a8 [label = 'Intermediate DBS result &\nNo Blood Plasma result\n n = 11']

a9 [label = 'Households \n n = 2564\nIndividuals \n n = 4433']

a11 [label = 'Households \n n = 1802\nIndividuals \n n = 2768']


a13 [label = 'Households \n n = 808\nIndividuals \n n = 1263']



data4 [label = 'Serology Analysis Data \n Baseline + Round 2 serology available\nNew positives | Baseline serology results', shape = folder, fillcolor = SpringGreen]
 


data6 [label = 'Leisuretime activity Data\nNew positives | Baseline serology results', shape = folder, fillcolor = YellowGreen]
data7 [label = 'Risk behavior Data\nNew positives | Baseline serology results', shape = folder, fillcolor = YellowGreen]

# edge definitions with the node IDs
data1-> a1 ;


a1-> data4[label='Round 2 DBS/Blood\npresent'];

a1 ->a6[label='Round 2 DBS/Blood missing'];

data4->a7
a7 -> a8[label='Intermediate\nresults excluded'];

a7 -> a9 ;
a9 -> data7[label='Risk behavior questions\ncompleted'];
a9 -> data6[xlabel='Leisuretime activity questions\ncompleted'];

data6->a13;
data7->a11;
{rank = same; a7; a8;}
{rank = same; a6; a1;}
                     

{rank = same; data6; data7}
}")
t
svg <- export_svg(t)
library(htmltools)
html_print(HTML(svg))

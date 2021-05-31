##### SEIR Diagram
library(groundhog)
groundhog.day <- as.Date("2021-01-01")
groundhog.library(DiagrammeR, groundhog.day)
groundhog.library(DiagrammeRsvg, groundhog.day)
groundhog.library(magrittr, groundhog.day)
groundhog.library(rsvg, groundhog.day)

SEIR_Digraph <- grViz("
  
  digraph boxes_and_circles {
    
    # Add node statements
    node [shape = box]
    S; E; R;
    Ip [label =<I<SUP>p</SUP>>];
    Ia [label =<I<SUP>&alpha;</SUP>>];
    Ik [label =<I<SUP>k</SUP>>];
    It1 [label =<I<SUP>t1</SUP>>];
    It2 [label =<I<SUP>t2</SUP>>];
    In [label =<I<SUP>n</SUP>>];
    
    # Add edge statements
    edge [color = black]
    S->E; 
    E->Ia [label = 'f']; 
    E->Ip [label = '1 - f'];
    Ia->R;
    Ip->It1
    Ip->Ik;
    Ip->In;
    It1->It2;
    It2->R;
    Ik->R;
    In->R;
  }


")

### Tried exporting it this way but it doesn't work, just cheated and used the Rstudio GUI instead
#grViz(SEIR_Digraph) %>% export_svg %>% charToRaw %>% rsvg_pdf("outputs/SEIR_Diagram.pdf")


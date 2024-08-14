import drilling_flow_cycle
import profile

cycle2 = [(0, 4),
 (0, 5),
 (2, 1),
 (2, 2),
 (1, 5),
 (0, 5),
 (2, 4),
 (2, 2),
 (1, 0),
 (0, 4),
 (0, 5),
 (2, 1),
 (2, 2),
 (1, 5),
 (0, 5),
 (2, 4),
 (2, 2),
 (1, 0),
 (0, 4)]

# profile.run("drilling_flow_cycle.drill_flow_cycle('dLQbccchhfo_122', cycle2)")   

cycle = [(2, 1), (5, 4), (4, 1), (7, 0), (6, 1)]

profile.run("drilling_flow_cycle.drill_flow_cycle('jLALPLQcbbegfihiihhrwhhaawr_122110011', cycle)")
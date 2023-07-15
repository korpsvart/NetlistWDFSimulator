function [thevImp] = adaptRJunction(Gprime, ZParams, refEdgeEndpoints, graphNodesOrder)

%Return the free parameter value needed to adapt a port of
%a generic R-type junction
%This is same as finding the equivalent Thevenin Impedance


%Inputs
%-Gprime: graph without the edge of the port to adapt (reference edge)
%-Zparams: already computed Z parameters
%-refEdgeEndpoints: endpoints (nodes) of the reference edge
%-graphNodesOrder: list of graph nodes, ordered as they are in the
%                  graph object


%Remove one of the endpoints of the reference edge
i1 = find(graphNodesOrder == refEdgeEndpoints(1));
i2 =find(graphNodesOrder == refEdgeEndpoints(2));

%Take them in this order so referenceRow < removeRow
%(no need to shift matrix rows in this way)
endpointToRemove = max([i1, i2]);
referenceEndpoint = min([i1, i2]);

%Get adjacency matrix
A = full(incidence(Gprime));

%We can't use the rmnode function because it will also delete incident
%edges, which is not what we want
%Instead we must remove the corresponding row. Hopefully matlab sorts
%the row according to node indexes, so we can exploit the natural index
%ordering...
A(endpointToRemove, :) = [];

%Get G conductance vector (and diag matrix)
%We've the Z vector already at our disposal so it's easy
G_vector = 1./ZParams;
G_m = diag(G_vector);

%Compute Y
Y = A*G_m*A'; %Admittance matrix
Imp = inv(Y); %Impedances matrix
%Adaptation
thevImp = Imp(referenceEndpoint, referenceEndpoint);



end


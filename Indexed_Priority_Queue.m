classdef Indexed_Priority_Queue < handle
    
    properties
        tree
        index
    end
    
    methods
        function ipq = Indexed_Priority_Queue(tree)
            % Constructor
            [~,order] = sort(tree(2,:)); % sort on taus
            ipq.tree = tree(:,order); % tree is a 2xM array with rxn # in row 1 and tau in row 2
            [~,ipq.index] = sort(ipq.tree(1,:)); % position of rxn i in tree
        end
        
        function [u,tau] = get_minimum(ipq)
            % Gets top node of priority queue, returns index of reaction u
            % and time tau
            u = ipq.tree(1,1);
            tau = ipq.tree(2,1);
        end
        
        function tau = get_rxn(ipq,rxn)
            % Returns time tau of rxn
            tau = get_value(ipq,ipq.index(rxn));
        end
        
        function update_rxn(ipq,rxn,new_tau)
            % Changes the tau in reaction rxn
            % Looks up position i in tree for rxn
            update(ipq,ipq.index(rxn),new_tau);
        end
    end
end

function swap(ipq,i,j)
    % i,j are nodes in the tree
    % Swap nodes
    temp = ipq.tree(:,i);
    ipq.tree(:,i) = ipq.tree(:,j);
    ipq.tree(:,j) = temp;
    % Update index
    temp = ipq.index(ipq.tree(1,i));
    ipq.index(ipq.tree(1,i)) = ipq.index(ipq.tree(1,j));
    ipq.index(ipq.tree(1,j)) = temp;
end

function update(ipq,i,new_val)
    % Changes value in node i, updates tree
    set_value(ipq,i,new_val);
    update_aux(ipq,i);
end

% RECURSION NOT ALLOWED IN MATLAB CODER!!!
% function update_aux(ipq,i)
%     % Rearrange tree in response to a changed node i
%     
%     parent = get_parent(ipq,i);
%     if has_children(ipq,i)
%         min_child = get_min_child(ipq,i);
%     end
%     
%     if get_value(ipq,i) < get_value(ipq,parent)
%         swap(ipq,i,parent);
%         update_aux(ipq,parent);
%     elseif has_children(ipq,i) && get_value(ipq,i) > get_value(ipq,min_child)
%         swap(ipq,i,min_child);
%         update_aux(ipq,min_child);
%     else
%         % done
%     end
% end

function update_aux(ipq,i)
    % Rearrange tree in response to a changed node i
    
    while 1 % assume node i is in wrong place
        
        parent = get_parent(ipq,i);
        if has_children(ipq,i)
            min_child = get_min_child(ipq,i);
        else
            min_child = 0;
        end
        
        if get_value(ipq,i) < get_value(ipq,parent)
            swap(ipq,i,parent);
            i = parent;
        elseif has_children(ipq,i) && get_value(ipq,i) > get_value(ipq,min_child)
            swap(ipq,i,min_child);
            i = min_child;
        else
            break % node i in right place
        end
        
    end
    

end

function set_value(ipq,i,val)
    % Set value tau of node i
    % Doesn't rebalance tree
    ipq.tree(2,i) = val;
end

function val = get_value(ipq,i)
    % Get value tau of node i
    val = ipq.tree(2,i);
end

function par = get_parent(ipq,i)
    % Gets position of parent of node i
    %   Returns own position if top node
    if i == 1
        par = 1;
    else
        par = floor(i/2);
    end
end

function hc = has_children(ipq,i)
    % Whether node i has children (1 or 2)
    if size(ipq.tree,2) > 2*i
        hc = 1;
    else
        hc = 0;
    end
end

function chs = get_children(ipq,i)
    % Gets positions of children (or child) of node i
    if size(ipq.tree,2) == 2*i % 1 child
        chs = 2*i;
    else
        chs = [2*i,2*i+1];
    end
end

function mch = get_min_child(ipq,i)
    % Gets position of child of i with minimum value (could be 1 child)
    chs = get_children(ipq,i);
    if length(chs) == 1 || get_value(ipq,chs(1)) < get_value(ipq,chs(2))
        mch = chs(1);
    else
        mch = chs(2);
    end
end

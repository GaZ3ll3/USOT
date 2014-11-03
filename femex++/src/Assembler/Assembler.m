classdef Assembler < handle
%DATABASE Hypothetical Matlab database API.
properties (Access = private)
  id_ % ID of the session instance.
end

methods
  function this = Assembler()
  %DATABASE Create a new database.
    this.id_ = Assembler_('new');
  end

  function delete(this)
  %DELETE Destructor.
    Assembler_('delete', this.id_);
  end
  
  function [F, DX, DY] = reference(this, nodes, qnodes)
  	[F, DX, DY] = Assembler_('reference', this.id_, nodes, qnodes);
  end
  
  function [I, J, V] = assema(this, pnodes, pelems, ref_fnk, weights, extern)
  	[I, J, V] = Assembler_('assema', this.id_, pnodes, pelems, ref_fnk, weights, extern);
  end
  
  function [I, J, V] = assems(this, pnodes, pelems, ref_gradx, ref_grady, weights, extern)
  	[I, J, V] = Assembler_('assems', this.id_, pnodes, pelems, ref_gradx, ref_grady, weights, extern);
  end
  
  function [L] = asseml(this, pnodes, qnodes, pelems, ref, weights, extern) 
    [L] = Assembler_('asseml', this.id_, pnodes, qnodes, pelems, ref, weights, extern);
  end
  
  function [C] = qnodes(this, pnodes, qnodes, pelems) 
  	[C] = Assembler_('qnodes', this.id_, pnodes, qnodes, pelems);
  end
  
  % Other methods...
end
end

classdef Boundary < handle
%DATABASE Hypothetical Matlab database API.
properties (Access = private)
  id_ % ID of the session instance.
end

methods
  function this = Boundary(type , edge)
  %DATABASE Create a new database.
    this.id_ = Boundary_('new', type, edge);
  end

  function delete(this)
  %DELETE Destructor.
    Boundary_('delete', this.id_);
  end
    
  
  function export(this)
      Boundary_('export', this.id_);
  end
  % Other methods...
end
end
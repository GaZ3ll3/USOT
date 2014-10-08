classdef Mesh < handle
%DATABASE Hypothetical Matlab database API.
properties (Access = private)
  id_ % ID of the session instance.
end

methods
  function this = Mesh(boundary)
  %DATABASE Create a new database.
    assert(isnumeric(boundary));
    this.id_ = Mesh_('new', boundary);
  end

  function delete(this)
  %DELETE Destructor.
    Mesh_('delete', this.id_);
  end
  
  function clear(this)
      Mesh_('clear', this.id_);
  end

  function report(this)
      Mesh_('report', this.id_);
  end

  % Other methods...
end
end

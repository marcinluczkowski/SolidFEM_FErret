using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SolidFEM.Classes
{
    public class PartClass
    {
        public List<Mesh> mesh { get; set; }
        public List<Material> material{ get; set; }

        public PartClass()
        {
            // empty
        }

        public PartClass(List<Mesh> _meshlist, List<Material> _material)
        {
            mesh = _meshlist;
            material = _material;
        }

    }

   
}

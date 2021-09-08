using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;

namespace SolidFEM.Classes
{
    public class Support
    {
        // -- properties --
        public Point3d Position { get; set; }
        public bool Tx { get; set; }
        public bool Ty { get; set; }
        public bool Tz { get; set; }

        // -- constructors --
        public Support()
        {
           // empty
        }

        public Support(Point3d _pos, bool _tx, bool _ty, bool _tz)
        {
            Position = _pos;
            Tx = _tx;
            Ty = _ty;
            Tz = _tz;
        }
        // -- methods --


    }
}

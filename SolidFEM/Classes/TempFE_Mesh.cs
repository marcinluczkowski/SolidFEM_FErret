using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;

namespace SolidFEM.Classes
{
    public class TempFE_Mesh
    {
        // -- parameters --
        public List<Mesh> MeshList { get; set; }
        public List<Node> MeshNodes { get; set; }
        public List<Element> MeshElements { get; set; }

        public List<double> MisesStress { get; private set; } // should not be able to modify after the analysis is done
        public List<double> dU { get; private set; }
        public List<double> dV { get; private set; }
        public List<double> dW { get; private set; }

        public Material Material { get; set; }

        // -- constructor --
        public TempFE_Mesh()
        {
            // empty
        }

        public TempFE_Mesh(List<Mesh> meshLst)
        {
            MeshList = meshLst;
        }

        public TempFE_Mesh(List<Mesh> meshLst, List<Node> nodes, List<Element> meshEls, List<double> mises, List<double> du, List<double> dv, List<double> dw, Material material)
        {
            MeshList = meshLst;
            MeshNodes = nodes;
            MeshElements = meshEls;
            MisesStress = mises;
            dU = du;
            dV = dv;
            dW = dw;
            Material = material;
        }

        // -- methods --
        /// <summary>
        /// Updates the displacements in the U-direction.
        /// </summary>
        /// <param name="newU"></param>
        public void UpdateUDisplacement(List<double> newU)
        {
            dU = newU;
        }
            
        public TempFE_Mesh DeepCopy()
        {
           return new TempFE_Mesh(this.MeshList, this.MeshNodes, this.MeshElements, this.MisesStress, this.dU, this.dV, this.dW, this.Material);
            
        }
    }
}

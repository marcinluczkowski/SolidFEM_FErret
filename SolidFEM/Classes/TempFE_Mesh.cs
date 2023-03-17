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
        public List<Mesh> MeshList { get; set; } // Remove
        public List<Node> MeshNodes { get; set; }
        public List<Element> MeshElements { get; set; }

        public List<double> Sigma_xx { get; set; }
        public List<double> Sigma_yy { get; private set; }
        public List<double> Sigma_zz { get; private set; }
        public List<double> Sigma_xy { get; private set; }
        public List<double> Sigma_xz { get; private set; }
        public List<double> Sigma_yz { get; private set; }
        public List<double> NodelMisesStresses { get; protected set; }
        public List<double> MisesStress { get; protected set; } // should not be able to modify after the analysis is done
        public List<double> dU { get; private set; }
        public List<double> dV { get; private set; }
        public List<double> dW { get; private set; }

        public MaterialOrto Material { get; set; }

        // -- constructor --
        public TempFE_Mesh()
        {
            // empty
        }

        /// <summary>
        /// Instantiate a FE Mesh from the list of meshes
        /// </summary>
        /// <param name="meshList">List of Rhino mesh.</param>
        public TempFE_Mesh(List<Mesh> meshList)
        {
            MeshList = meshList;
        }

        public TempFE_Mesh(List<Mesh> meshLst, List<Node> nodes, List<Element> meshEls, List<double> mises, List<double> nodalMises, List<double> du, List<double> dv, List<double> dw, MaterialOrto material,
            List<double> s_xx, List<double> s_yy, List<double> s_zz, List<double> s_xy, List<double> s_xz, List<double> s_yz)
        {
            MeshList = meshLst;
            MeshNodes = nodes;
            MeshElements = meshEls;
            MisesStress = mises;
            NodelMisesStresses = nodalMises;
            Sigma_xx = s_xx; Sigma_yy = s_yy; Sigma_zz = s_zz;
            Sigma_xy = s_xy; Sigma_xz = s_xz; Sigma_yz = s_yz;
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
           return new TempFE_Mesh(this.MeshList, this.MeshNodes, this.MeshElements, this.MisesStress, this.NodelMisesStresses ,this.dU, this.dV, this.dW, this.Material, 
               this.Sigma_xx, this.Sigma_yy, this.Sigma_zz, this.Sigma_xy, this.Sigma_xz, this.Sigma_yz);
            
        }
    }
}

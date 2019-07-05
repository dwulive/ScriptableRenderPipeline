using UnityEditor.ShaderGraph;

namespace UnityEditor.Rendering.HDPipeline
{
    static class CreateDecalShaderGraph
    {
        [MenuItem("Assets/Create/Shader/HDRP/Decal Graph", false, 208)]
        public static void CreateMaterialGraph()
        {
            GraphUtil.CreateNewGraph(new DecalMasterNode());
        }
    }
}

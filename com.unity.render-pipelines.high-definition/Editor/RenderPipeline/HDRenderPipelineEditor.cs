using UnityEngine.Rendering.HDPipeline;

namespace UnityEditor.Rendering.HDPipeline
{
    [CustomEditor(typeof(HDRenderPipelineAsset))]
    [CanEditMultipleObjects]
    public sealed class HDRenderPipelineEditor : Editor
    {
        SerializedHDRenderPipelineAsset m_SerializedHDRenderPipeline;

        internal bool showInspector = true;

        void OnEnable()
        {
#if QUALITY_SETTINGS_GET_RENDER_PIPELINE_AT_AVAILABLE
            showInspector = false;
#endif
            m_SerializedHDRenderPipeline = new SerializedHDRenderPipelineAsset(serializedObject);
        }

        public override void OnInspectorGUI()
        {
            if (!showInspector)
                return;

            var serialized = m_SerializedHDRenderPipeline;

            serialized.Update();

            HDRenderPipelineUI.Inspector.Draw(serialized, this);

            serialized.Apply();
        }
    }
}

using UnityEngine.Rendering.HDPipeline;
using UnityEngine.Rendering;

namespace UnityEditor.Rendering.HDPipeline
{
    [CanEditMultipleObjects]
    [VolumeComponentEditor(typeof(VolumetricLightingController))]
    public class VolumetricLightingControllerEditor : VolumeComponentEditor
    {
        public override void OnInspectorGUI()
        {
            base.OnInspectorGUI();

            if (!(GraphicsSettings.renderPipelineAsset as HDRenderPipelineAsset)
                ?.currentPlatformRenderPipelineSettings.supportVolumetrics ?? false)
            {
                EditorGUILayout.Space();
                EditorGUILayout.HelpBox("The current HDRP Asset does not support Volumetrics.", MessageType.Error, wide: true);
            }
        }
    }
}


#include "io_environment.h"

#include "sph_system.h"

#ifdef _WIN32
#define PATHSEP '\\'
#else
#define PATHSEP '/'
#endif

namespace SPH
{
//=============================================================================================//
IOEnvironment::IOEnvironment(SPHSystem &sph_system, const std::string &workPath)
    : sph_system_(sph_system),
      input_folder_(workPath + PATHSEP + "input"), output_folder_(workPath + PATHSEP + "output"),
      restart_folder_(workPath + PATHSEP + "restart"), reload_folder_(workPath + PATHSEP + "reload")
{
    if (!fs::exists(input_folder_))
    {
        fs::create_directory(input_folder_);
    }

    if (!fs::exists(output_folder_))
    {
        fs::create_directory(output_folder_);
    }
    else
    {
        std::string backup_name = output_folder_ + "_backup";
        fs::path output_dir = output_folder_;
        //fs::path backup_path = output_dir.parent_path() / backup_name;
        //if (fs::exists(backup_path))
        //{
        //    fs::remove_all(backup_path);
        //}
        //fs::rename(output_dir, backup_path); // cause crash in Windows
        //std::cout << "Moved existing output to: " << backup_path << std::endl;
        fs::create_directory(output_folder_);
    }

    if (!fs::exists(restart_folder_))
    {
        fs::create_directory(restart_folder_);
    }

    if (!fs::exists(reload_folder_))
    {
        fs::create_directory(reload_folder_);
    }

    sph_system.io_environment_ = this;
}
//=================================================================================================//
void IOEnvironment::resetForRestart()
{
    fs::remove_all(restart_folder_);
    fs::create_directory(restart_folder_);
}
//=================================================================================================//
void IOEnvironment::reinitializeReloadFolder()
{
    fs::remove_all(reload_folder_);
    fs::create_directory(reload_folder_);
}
//=================================================================================================//
void IOEnvironment::appendOutputFolder(const std::string &append_name)
{
    output_folder_ += "_" + append_name;
    if (!fs::exists(output_folder_))
    {
        fs::create_directory(output_folder_);
    }
    else
    {
        fs::remove_all(output_folder_);
        fs::create_directory(output_folder_);
    }
}
//=================================================================================================//
void IOEnvironment::resetOutputFolder(const std::string &new_name)
{
    output_folder_ = new_name;
    if (!fs::exists(output_folder_))
    {
        fs::create_directory(output_folder_);
    }
    else
    {
        fs::remove_all(output_folder_);
        fs::create_directory(output_folder_);
    }
}
//=============================================================================================//
ParameterizationIO *IOEnvironment::defineParameterizationIO()
{
    return parameterization_io_ptr_keeper_.createPtr<ParameterizationIO>(input_folder_);
}
//=================================================================================================//
} // namespace SPH

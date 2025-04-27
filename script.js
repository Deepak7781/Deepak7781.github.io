const projects = [
    {
        title: "NDI Controller for a fighter aircraft",
        description: "Creating a NDI Controller for F-16 Fighter aircraft",
        link: "https://github.com/Kavin-Senthilkumar/NDI_controller_for_a_fighter_aircraft"
    },
    {
        title: "Blade Element Theory",
        description: "A simple simulation of Blade element theory using MATLAB",
        link: "https://github.com/Deepak7781/Blade_Element_Theory"
    }

];

function displayProjects(){
    const projectList = document.getElementById("project-list");
    projectList.innerHTML = "";

    projects.forEach(project => {
        const projectCard = document.createElement("div");
        projectCard.classList.add("project-card");

        projectCard.innerHTML = `
            <h3>${project.title}</h3>
            <p>${project.description}</p>
            <a href="${project.link}" target="_blank">View on GitHub</a>   
         `;

         projectList.appendChild(projectCard);
    });
}

window.onload = displayProjects;
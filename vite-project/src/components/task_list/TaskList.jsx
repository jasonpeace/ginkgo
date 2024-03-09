
import { useState, useContext, createContext } from 'react';

const ActionContext = createContext();

const TaskList = () => {
    const [tasks, setTasks] = useState(["Task1", "Task2", "Task3"]);


    const handleDelete = (item) => {
        let new_tasks = tasks.filter((task)=> task != item);
        setTasks(new_tasks)

    };

    const actions = {
            handleDelete:handleDelete
        }

    return(
        <ActionContext.Provider value={actions}>
        <div >
            <AddTaskComp tasks={tasks} handleAdd={setTasks}/>
            <TaskListComp tasks={tasks} />
        </div>
        </ActionContext.Provider>

    );
};



const AddTaskComp = ({tasks, handleAdd}) =>{

    const handleClick = (event) =>{
        let task = document.getElementById("addTaskField").value;
        if (!tasks){
            tasks = [];
        }
        let items = [...tasks, task]
        handleAdd(items)
    };

    return(
        <div>
                <input type="text" id="addTaskField" name="assTaskField" />
                <button onClick={handleClick}>Add Task </button>
        </div>
    );
};

const TaskListComp = ({tasks}) =>{

    return(
        <ul>
          {tasks.map((task)=>{
            return <TaskComp key={task} taskLabel={task} />;
            })
          }
        </ul>
    );
};


const TaskComp = ({taskLabel}) =>{
    const actionContext = useContext(ActionContext);
    const handleDelete = actionContext.handleDelete;

    const handleClick = () => {
        handleDelete(taskLabel)
    }

    return(
        <li>
            <label>{taskLabel}</label>
            <button onClick={handleClick}>remove</button>
        </li>

    );
};


export default TaskList;

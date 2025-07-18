This short guide has been written to facilitate contributions to the present Wiki.
We encourage Hibridon users to directly add/edit Wiki pages, or, report any missing or incorrect documentation using the [Issues board](https://github.com/hibridon/hibridon/issues) of the Hibridon repository.


1. [How do I create a page?](#how-do-i-create-a-page)
2. [How do I edit a page?](#how-do-i-edit-a-page)
3. [How do I insert a new line?](#how-do-i-insert-a-new-line)
4. [How do I create a link?](#how-do-i-create-a-link)
5. [How do I rename a page?](#how-do-i-rename-a-page)
6. [How can I add pictures to pages?](#how-can-i-add-pictures-to-pages)
7. [How can I write equations?](#how-can-i-write-equations)
8. [Are there any tools to make editing easier/faster?](#are-there-any-tools-to-make-editing-easierfaster)
9. [How can I update the sidebar menu?](#how-can-i-update-the-sidebar-menu)

## How do I create a page?

To create a new page, just click on the "New Page" green button locate at the top right of the page:
<p align="center"><img width="168" alt="Screenshot 2022-05-24 at 19 36 57" src="https://user-images.githubusercontent.com/20931653/170148861-6b09a3da-d167-468a-b277-b042f5a95286.png"></p>

You can then choose a title for this new page:
<p align="center"><img width="500" alt="Screenshot 2022-05-24 at 19 38 19" src="https://user-images.githubusercontent.com/20931653/170148998-2946996a-1718-460f-af41-d805f45ccd51.png"></p>


ℹ️ Try to choose something short but meaningful as the title of the page is used to create the associated markdown file that can then be referred to by [creating a link](#how-do-i-create-a-link).


⚠️ If you want your new page to appear as a subpage of another page in the sidebar menu, its title must start by the title of the parent page, followed by three spaces, and then followed by the title of the subpage:
  ```markdown
  Title of the parent page   Title of the subpage
  ```

## How do I edit a page?

Navigate to the Wiki page you want to edit and then click on the "Edit" button:
<p align="center"><img width="168" alt="Screenshot 2022-05-24 at 19 36 57" src="https://user-images.githubusercontent.com/20931653/170148861-6b09a3da-d167-468a-b277-b042f5a95286.png"></p>

You can then edit the page using Markdown syntax and/or HTML tags:
<p align="center"><img width="500" alt="Screenshot 2022-05-24 at 19 51 20" src="https://user-images.githubusercontent.com/20931653/170150231-de42be58-9cfc-4f10-abfc-110bbafca1c8.png"></p>


Before saving your modifications by clicking on the "Save Page" green button at the bottom of the editing page, you can preview the page by clicking on "Preview":
<p align="center"><img width="500" alt="Screenshot 2022-05-24 at 19 51 22" src="https://user-images.githubusercontent.com/20931653/170150262-00cfdf95-cf1e-4a17-b87e-f5f2b78fa400.png"></p>


## How do I insert a new line?

Two new lines in the edit mode causes a paragraph break in the final rendering:
  Markdown code:
  ```markdown
  This is a 
  
  paragraph break
  ```
  Result:
  >This is a
  >
  >paragraph break


## How do I create a link?

To create a link from a Wiki page to another within the edit mode, you can use the Markdown syntax:
  Markdown code:
  ```markdown
  [Text to print](url-of-the-page)
  ```
  Result:
  >[Test to print](url-of-the-page)

The url of a Wiki page is simply the title of the page, with all spaces replaced by dashes -, and all special characters (e.g. ! or ?) ignored.


ℹ️ It is also possible to use the HTML tag \<a> \</a> to create links.

## How do I rename a page?

While editing a Wiki page, you can change the title at the top of the editing page.

⚠️ Renaming the title of a page will also change its url. All links to this page will then need to be updated with a the url. This operation can be done automatically by cloning the Wiki repository locally and using grep  and sed or a code editor to search and replace the links url by the correct one. If your are not sure, please use the [Issues board](https://github.com/hibridon/hibridon/issues) of the Hibridon directory to ask someone to rename the page.

## How can I add pictures to pages?

While in edit mode, you can directly add a picture dragging a file (png, jpg, gif, ...) from your local file explorer to the text area of the edit page.
Your file will automatically be uploaded with a specific url and the something similar to the following will appear in the edit zone:
```markdown
![Image Name](https://user-images.githubusercontent.com/20931653/170152147-2f0da3a1-78d4-43ff-ae4c-4115cb5sdd79.png)
```

It is not possible to resize pictures added this way. If you want more control on how the image will appear you can copy the image url and use the \<img> HTML tag:
```html
<img height="200px" width="150px" src="https://user-images.githubusercontent.com/20931653/170152147-2f0da3a1-78d4-43ff-ae4c-4115cb5sdd79.png">
```

## How can I write equations?

GitHub just (half of may 2022) added support of LaTeX in the rendering of markdown files present in the repositories. Unfortunately, this features is not yet enabled on the Wiki pages.
This should be corrected soon.

Right now, the only solution is to include images of the equations. This can be easily done using tools such as `LaTeXiT`.

## Are there any tools to make editing easier/faster?

Using Git, this Wiki repository can be cloned locally:
<p align="center"><img width="378" alt="Screenshot 2022-05-24 at 20 20 48" src="https://user-images.githubusercontent.com/20931653/170152691-c5385f9c-0d95-4f41-b209-d96a0206bfd5.png"></p>

Once cloned on your computer, you can manually edit any markdown file, commit the modifications and then push them using Git.

## How can I update the sidebar menu?

The sidebar menu is automatically generated and it updates itself at each modification of any of the pages using the [github-wiki-sidebar](https://www.npmjs.com/package/github-wiki-sidebar) tool. It may take up to a few minutes for the updated sidebar menu to appear.

⚠️ Please do not edit the sidebar (`_sidebar.md`) manually. 
